import numpy as np
from scipy.sparse import coo_matrix
from tqdm import tqdm


def make_kmer_sparse_matrix(sigobjs):
    print("Making sparse k-mer matrix")

    column_indices = []
    row_indices = []
    sparse_data = []

    hashval_to_index = {}

    len(hashval_to_index)

    for i, sig in tqdm(enumerate(sigobjs)):

        # Get hashval position
        if i == 0:
            hashval_to_index.update(
                dict({hashval: j for j, hashval in enumerate(sig.minhash.hashes)})
            )
        else:
            # Next column position is the length (since indices are 1-based)
            j = len(hashval_to_index)
            for hashval in sig.minhash.hashes:
                if hashval not in hashval_to_index:
                    hashval_to_index[hashval] = j
                    j += 1
        sig_col = [hashval_to_index[hashval] for hashval in sig.minhash.hashes]
        sig_row = [i] * len(sig_col)
        sig_data = [abundance for abundance in sig.minhash.hashes.values()]

        column_indices.extend(sig_col)
        row_indices.extend(sig_row)
        sparse_data.extend(sig_data)

    kmer_coo_matrix = coo_matrix((sparse_data, (row_indices, column_indices)))
    return kmer_coo_matrix, hashval_to_index


def sparse_var(a, axis=None):
    """Variance of sparse matrix a
    var = mean(a**2) - mean(a)**2
    """
    a_squared = a.copy()
    a_squared.data **= 2
    return a_squared.mean(axis) - np.square(a.mean(axis))


def sparse_std(a, axis=None):
    """Standard deviation of sparse matrix a
    std = sqrt(var(a))
    """
    return np.sqrt(vars(a, axis))
