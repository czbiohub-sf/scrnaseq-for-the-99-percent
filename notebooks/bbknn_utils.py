from typing import Optional, Union


import numpy as np
from scipy.spatial import cKDTree
import bbknn
from scipy.sparse import csr_matrix
import scanpy as sc
from numpy.testing import assert_array_equal, assert_array_compare
import operator


import numpy as np
from anndata import AnnData
from sklearn.utils import check_random_state, check_array

from scanpy.tools._utils import get_init_pos_from_paga#, _choose_representation
from scanpy import logging as logg
from scanpy._settings import settings
from scanpy._compat import Literal
from scanpy._utils import AnyRandom, NeighborsView


# Lots of this was stolen from https://github.com/theislab/scanpy/blob/master/scanpy/tools/_umap.py
_InitPos = Literal['paga', 'spectral', 'random']


def make_graph_from_batch_corrected_distances(distances, batch_list, neighbors_within_batch, approx, metric, use_faiss, n_trees):
    '''
    Identify the KNN structure to be used in graph construction. All input as in ``bbknn.bbknn()``
    and ``bbknn.bbknn_pca_matrix()``. Returns a tuple of distances and indices of neighbours for
    each cell.
    '''
    #get a list of all our batches
    batches = np.unique(batch_list)
    #in case we're gonna be faissing, turn the data to float32
    if metric=='euclidean' and not approx and 'faiss' in sys.modules and use_faiss:
        pca = pca.astype('float32')
    #create the output matrices, with the indices as integers and distances as floats
    knn_distances = np.zeros((distances.shape[0],neighbors_within_batch*len(batches)))
    knn_indices = np.copy(knn_distances).astype(int)
    #find the knns using faiss/cKDTree/KDTree/annoy
    #need to compare each batch against each batch (including itself)
    for to_ind in range(len(batches)):
        #this is the batch that will be used as the neighbour pool
        #create a boolean mask identifying the cells within this batch
        #and then get the corresponding row numbers for later use
        batch_to = batches[to_ind]
        mask_to = batch_list == batch_to
        ind_to = np.arange(len(batch_list))[mask_to]
        #create the faiss/cKDTree/KDTree/annoy, depending on approx/metric
        ckd = bbknn.create_tree(data=distances[mask_to,:],approx=approx,metric=metric,
                          use_faiss=use_faiss,n_trees=n_trees)
        for from_ind in range(len(batches)):
            #this is the batch that will have its neighbours identified
            #repeat the mask/row number getting
            batch_from = batches[from_ind]
            mask_from = batch_list == batch_from
            ind_from = np.arange(len(batch_list))[mask_from]
            #fish the neighbours out, getting a (distances, indices) tuple back
            ckdout = bbknn.query_tree(data=distances[mask_from,:],ckd=ckd,
                                neighbors_within_batch=neighbors_within_batch,
                                approx=approx,metric=metric,use_faiss=use_faiss)
            #the identified indices are relative to the subsetted PCA matrix
            #so we need to convert it back to the original row numbers
            for i in range(ckdout[1].shape[0]):
                for j in range(ckdout[1].shape[1]):
                    ckdout[1][i,j] = ind_to[ckdout[1][i,j]]
            #save the results within the appropriate rows and columns of the structures
            col_range = np.arange(to_ind*neighbors_within_batch, (to_ind+1)*neighbors_within_batch)
            knn_indices[ind_from[:,None],col_range[None,:]] = ckdout[1]
            knn_distances[ind_from[:,None],col_range[None,:]] = ckdout[0]
    return knn_distances, knn_indices


def bbknn_distance_matrix(distances, batch_list, neighbors_within_batch=3, trim=None, 
          approx=True, n_trees=10, use_faiss=True, metric='angular', 
          set_op_mix_ratio=1, local_connectivity=1):
    '''
    Scanpy-independent BBKNN variant that runs on a PCA matrix and list of per-cell batch assignments instead of
    an AnnData object. Non-data-entry arguments behave the same way as ``bbknn.bbknn()``.
    Returns a ``(distances, connectivities)`` tuple, like what would have been stored in the AnnData object.
    The connectivities are the actual neighbourhood graph.
    
    Input
    -----
    pca : ``numpy.array``
        PCA (or other dimensionality reduction) coordinates for each cell, with cells as rows.
    batch_list : ``numpy.array`` or ``list``
        A list of batch assignments for each cell.
    '''
    #more basic sanity checks/processing
    #do we have the same number of cells in pca and batch_list?
    if distances.shape[0] != len(batch_list):
        raise ValueError("Different cell counts indicated by `distances.shape[0]` and `len(batch_list)`.")
    #convert batch_list to np.array of strings for ease of mask making later
    batch_list = np.asarray([str(i) for i in batch_list])
    #metric sanity checks (duplicating the ones in bbknn(), but without scanpy logging)
    if approx and metric not in ['angular', 'euclidean', 'manhattan', 'hamming']:
        print('unrecognised metric for type of neighbor calculation, switching to angular')
        metric = 'angular'
    elif not approx and not (metric=='euclidean' or isinstance(metric,DistanceMetric) or metric in KDTree.valid_metrics):
        print('unrecognised metric for type of neighbor calculation, switching to euclidean')
        metric = 'euclidean'
    #obtain the batch balanced KNN graph
    knn_distances, knn_indices = make_graph_from_batch_corrected_distances(
        distances, 
        batch_list=batch_list, 
        n_trees=n_trees,
        approx=approx,
        metric=metric,
        use_faiss=use_faiss,
        neighbors_within_batch=neighbors_within_batch)
    #sort the neighbours so that they're actually in order from closest to furthest
    newidx = np.argsort(knn_distances,axis=1)
    knn_indices = knn_indices[np.arange(np.shape(knn_indices)[0])[:,np.newaxis],newidx]
    knn_distances = knn_distances[np.arange(np.shape(knn_distances)[0])[:,np.newaxis],newidx] 
    #this part of the processing is akin to scanpy.api.neighbors()
    dist, cnts = bbknn.compute_connectivities_umap(knn_indices, knn_distances, knn_indices.shape[0], 
                                             knn_indices.shape[1], set_op_mix_ratio=set_op_mix_ratio,
                                             local_connectivity=local_connectivity)
    #trimming. compute default range if absent
    if trim is None:
        trim = 10 * knn_distances.shape[1]
    #skip trimming if set to 0, otherwise trim
    if trim > 0:
        cnts = bbknn.trimming(cnts=cnts,trim=trim)
    return (dist, cnts)


def assign_neighbors(ad, neighbors_key, knn_distances, knn_indices, set_use_rep=True):
    """Add bbknn-corrected neighbors to specific keybor key"""
    ad.uns[neighbors_key] = {}
    #we'll have a zero distance for our cell of origin, and nonzero for every other neighbour computed
    ad.uns[neighbors_key]['params'] = {
        'n_neighbors': len(knn_distances[0,:].data)+1, 
        'method': 'umap', 
        # Need this to force UMAP to use the raw data as the representation
         'use_rep': "X"
    }
    distances_key = f'{neighbors_key}__distances'
    connectivities_key = f'{neighbors_key}__connectivities'
    ad.obsp[connectivities_key] = csr_matrix(knn_indices)
    ad.obsp[distances_key] = knn_distances
    ad.uns[neighbors_key]['distances_key'] = distances_key
    ad.uns[neighbors_key]['connectivities_key'] = connectivities_key
#     ad.uns[neighbors_key]['distances'] = knn_distances
#     ad.uns[neighbors_key]['connectivities'] = csr_matrix(knn_indices)
    ad.uns[neighbors_key]['params'] = {}
    ad.uns[neighbors_key]['params']['metric'] = 'precomputed'
    ad.uns[neighbors_key]['params']['method'] = 'bbknn'
    if set_use_rep:
        ad.uns[neighbors_key]['params']['use_rep'] = neighbors_key

    return ad


def bbknn_similarity_matrix_and_assign_adata(
    adata,
    obsp_key,
    color=['narrow_group', 'species', 'PTPRC', 'SFTPC', 'n_counts', 'n_genes'],
    COUNTS_BASED_UMAP_COORDS=None,
    neighbors_within_batch=15,
    set_use_rep=True,
    batch_key='species'
    **kwargs,
):
    """
    COUNTS_BASED_UMAP_COORDS : array
        Used for a sanity check to assert that the new UMAP doesn't match the original one
    """
    print(adata)
    batch_list = adata.obs[batch_key].tolist()
    print(f"len(batch_list): {len(batch_list)}")
#     import pdb; pdb.set_trace()

    # Get similarity; stored at pairwise location
    data = adata.obsp[obsp_key]
    knn_distances, knn_indices = bbknn_distance_matrix(
        distances=data, batch_list=batch_list, neighbors_within_batch=neighbors_within_batch)
    
#     import pdb; pdb.set_trace()
    adata = assign_neighbors(adata, obsp_key, knn_distances, knn_indices, set_use_rep=set_use_rep)

    sc.tl.umap(adata, neighbors_key=obsp_key, **kwargs)
#     umap_precomputed(adata, neighbors_key=neighbors_key, **kwargs)
    
    if COUNTS_BASED_UMAP_COORDS is not None:
        assert_array_compare(operator.__ne__, COUNTS_BASED_UMAP_COORDS, adata.obsm['X_umap'])
 
    sc.pl.umap(adata, neighbors_key=obsp_key, color=color, ncols=2)



# def _choose_representation(adata, use_rep=None, n_pcs=None, silent=False):
#     verbosity = settings.verbosity
#     if silent and settings.verbosity > 1:
#         settings.verbosity = 1
#     if use_rep is None and n_pcs == 0:  # backwards compat for specifying `.X`
#         logg.warning('use_rep=None and n_pcs=0')
#         use_rep = 'X'
#     if use_rep is None:
#         logg.warning('use_rep=None')
#         if adata.n_vars > settings.N_PCS:
#             logg.warning('adata.n_vars > settings.N_PCS')
#             if 'X_pca' in adata.obsm.keys():
#                 if n_pcs is not None and n_pcs > adata.obsm['X_pca'].shape[1]:
#                     raise ValueError(
#                         '`X_pca` does not have enough PCs. Rerun `sc.pp.pca` with adjusted `n_comps`.')
#                 X = adata.obsm['X_pca'][:, :n_pcs]
#                 logg.info(f'    using \'X_pca\' with n_pcs = {X.shape[1]}')
#             else:
#                 logg.warning(
#                     f'You’re trying to run this on {adata.n_vars} dimensions of `.X`, '
#                     'if you really want this, set `use_rep=\'X\'`.\n         '
#                     'Falling back to preprocessing with `sc.pp.pca` and default params.'
#                 )
#                 X = pca(adata.X)
#                 adata.obsm['X_pca'] = X[:, :n_pcs]
#         else:
#             logg.info('    using data matrix X directly')
#             X = adata.X
#     else:
#         if use_rep in adata.obsm.keys():
#             X = adata.obsm[use_rep]
#             if use_rep == 'X_pca' and n_pcs is not None:
#                 X = adata.obsm[use_rep][:, :n_pcs]
#         elif use_rep == 'X':
#             X = adata.X
#         else:
#             raise ValueError(
#                 'Did not find {} in `.obsm.keys()`. '
#                 'You need to compute it first.'.format(use_rep))
#     settings.verbosity = verbosity  # resetting verbosity
#     return X
    

# def umap_precomputed(
#     adata: AnnData,
#     min_dist: float = 0.5,
#     spread: float = 1.0,
#     n_components: int = 2,
#     maxiter: Optional[int] = None,
#     alpha: float = 1.0,
#     gamma: float = 1.0,
#     negative_sample_rate: int = 5,
#     init_pos: Union[_InitPos, np.ndarray, None] = 'spectral',
#     random_state: AnyRandom = 0,
#     a: Optional[float] = None,
#     b: Optional[float] = None,
#     copy: bool = False,
#     method: Literal['umap', 'rapids'] = 'umap',
#     neighbors_key: Optional[str] = None,
#     COUNTS_BASED_UMAP_COORDS=None,
# ):
#     adata = adata.copy() if copy else adata

#     if neighbors_key is None:
#         neighbors_key = 'neighbors'

#     if neighbors_key not in adata.uns:
#         raise ValueError(
#             f'Did not find .uns["{neighbors_key}"]. Run `sc.pp.neighbors` first.')
#     start = logg.info('computing UMAP')

#     neighbors = NeighborsView(adata, neighbors_key)

#     if ('params' not in neighbors
#         or neighbors['params']['method'] != 'umap'):
#         logg.warning(f'.obsp["{neighbors["connectivities_key"]}"] have not been computed using umap')
#     from umap.umap_ import find_ab_params, simplicial_set_embedding
#     if a is None or b is None:
#         a, b = find_ab_params(spread, min_dist)
#     else:
#         a = a
#         b = b
#     adata.uns['umap'] = {'params':{'a': a, 'b': b}}
#     if isinstance(init_pos, str) and init_pos in adata.obsm.keys():
#         init_coords = adata.obsm[init_pos]
#     elif isinstance(init_pos, str) and init_pos == 'paga':
#         init_coords = get_init_pos_from_paga(adata, random_state=random_state, neighbors_key=neighbors_key)
#     else:
#         init_coords = init_pos  # Let umap handle it
#     if hasattr(init_coords, "dtype"):
#         init_coords = check_array(init_coords, dtype=np.float32, accept_sparse=False)

#     if random_state != 0:
#         adata.uns['umap']['params']['random_state'] = random_state
#     random_state = check_random_state(random_state)

#     neigh_params = neighbors['params']
#     X = _choose_representation(
#         adata, neigh_params.get('use_rep', None), neigh_params.get('n_pcs', None), silent=True)
    
#     # ---- debugger ---- #
# #     import pdb; pdb.set_trace()
#     # ---- debugger ---- #


#     if method == 'umap':
#         # the data matrix X is really only used for determining the number of connected components
#         # for the init condition in the UMAP embedding
#         n_epochs = 0 if maxiter is None else maxiter
#         X_umap = simplicial_set_embedding(
#             X,
#             neighbors['connectivities'].tocoo(),
#             n_components,
#             alpha,
#             a,
#             b,
#             gamma,
#             negative_sample_rate,
#             n_epochs,
#             init_coords,
#             random_state,
#             neigh_params.get('metric', 'euclidean'),
#             neigh_params.get('metric_kwds', {}),
#             verbose=settings.verbosity > 3,
#         )
#     elif method == 'rapids':
#         metric = neigh_params.get('metric', 'euclidean')
#         if metric != 'euclidean':
#             raise ValueError(
#                 f'`sc.pp.neighbors` was called with `metric` {metric!r}, '
#                 "but umap `method` 'rapids' only supports the 'euclidean' metric."
#             )
#         from cuml import UMAP
#         n_neighbors = neighbors['params']['n_neighbors']
#         n_epochs = 500 if maxiter is None else maxiter # 0 is not a valid value for rapids, unlike original umap
#         X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
#         umap = UMAP(
#             n_neighbors=n_neighbors,
#             n_components=n_components,
#             n_epochs=n_epochs,
#             learning_rate=alpha,
#             init=init_pos,
#             min_dist=min_dist,
#             spread=spread,
#             negative_sample_rate=negative_sample_rate,
#             a=a,
#             b=b,
#             verbose=settings.verbosity > 3,
#         )
#         X_umap = umap.fit_transform(X_contiguous)
#     adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
#     logg.info(
#         '    finished',
#         time=start,
#         deep=(
#             'added\n'
#             "    'X_umap', UMAP coordinates (adata.obsm)"
#         ),
#     )
#     return adata if copy else None
