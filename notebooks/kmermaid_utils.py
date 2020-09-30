import glob
import os
import sys

from joblib import Parallel, delayed
import pandas as pd
import sourmash
from sourmash.logging import error, notify
from sourmash import sourmash_args

from tqdm import tqdm



SKETCH_INFO_PATTERN = '(?P<mol_or_alpha>molecule|alphabet)-(?P<alphabet>\w+)__ksize-(?P<ksize>\d+)__(?P<sketch_style>num_hashes|scaled)-(?P<sketch_value>\d+)'

SKETCH_PARAMS = ['alphabet', 'ksize', 'sketch_style', 'sketch_value']

CELL_PARAMS = ['channel', 'cell_barcode']
SKETCH_CELL_PARAMS = SKETCH_PARAMS + CELL_PARAMS

def make_sketch_id(alpha, ksize, style, value):
    sketch_id = f'alphabet-{alpha}__ksize-{ksize}__{style}-{value}'
    return sketch_id


def get_sig_file_df(sig_folder, verbose=False, aligned_unaligned_merged=False):
    df = pd.Series(glob.glob(f'{sig_folder}/*.sig'), name='fullpath').to_frame()
    df['basename'] = df['fullpath'].map(os.path.basename)
    if verbose:
        print('df.shape:', df.shape)

    if aligned_unaligned_merged:
        cell_info = pd.Series(df['basename'].str.split('.').str[0], name='cell_id').to_frame()
        sig_info = pd.DataFrame(index=df.index)
    else:
        sig_info = df['basename'].str.extractall(SKETCH_INFO_PATTERN)
        sig_info = sig_info.droplevel(-1)
        sig_info['ksize'] = sig_info['ksize'].astype(int)
        sig_info['sketch_value'] = sig_info['sketch_value'].astype(int)
        cell_info = df['basename'].str.split('__', expand=True)
        cols_to_drop = cell_info.columns.intersection([3,4,5,6,7])
        cell_info = cell_info.drop(cols_to_drop, axis=1)
        cell_info = cell_info.rename(columns={
            0: 'channel', 1: 'alignment_status', 2: 'cell_barcode'})
        cell_info['cell_id'] = cell_info.apply(
            lambda x: "{channel}__{cell_barcode}".format(**x), axis=1)

    
    sigs_with_metadata = pd.concat([df, sig_info, cell_info], axis=1)
    
    # don't worrry about sorting for now
#     sort_cols = sigs_with_metadata.columns.intersection(SKETCH_CELL_PARAMS)
#     if len(sort_cols) > 0:
#         sigs_with_metadata = sigs_with_metadata.sort_values(sort_cols)
    if verbose:
        print('sigs_with_metadata.shape:', sigs_with_metadata.shape)
    return sigs_with_metadata

## -- Sourmash signature stuff --- ###


# from sourmash.sig import _check_abundance_compatibility

def merge(filenames, ksize, moltype, flatten=False):
    """
    merge one or more signatures.
    
    Adapted from 'sourmash sig merge' command
    """

    first_sig = None
    mh = None
    total_loaded = 0


    for sigfile in filenames:
        this_n = 0
        for sigobj in sourmash.load_signatures(sigfile,
                                                        ksize=ksize,
                                                        select_moltype=moltype):

            # first signature? initialize a bunch of stuff
            if first_sig is None:
                first_sig = sigobj
                mh = first_sig.minhash.copy_and_clear()

            try:
                sigobj_mh = sigobj.minhash

                mh.merge(sigobj_mh)
            except:
                error("ERROR when merging signature '{}' ({}) from file {}",
                      sigobj.name(), sigobj.md5sum()[:8], sigfile)
                raise

            this_n += 1
            total_loaded += 1

    merged_sigobj = sourmash.SourmashSignature(mh)
    
    return merged_sigobj



def merge_signatures(params, df, outdir_base):
    alpha, ksize, style, value, channel, barcode = params
    
    sketch_id = f'alphabet-{alpha}__ksize-{ksize}__{style}-{value}'
    outdir = f"{outdir_base}/{sketch_id}"
    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except FileExistsError:
            # Already created by another thread
            pass
    
    cell_id = f'{channel}__{barcode}'
    outsig = f'{outdir}/{cell_id}.sig'
    insigs = ' '.join(df['fullpath'])
    
    try:
        merged_sigobj = merge(df['fullpath'].values, ksize=ksize, moltype=alpha)
    except ValueError:
        print(f'cell: {cell_id} was not able to be merged\n')
        return cell_id, False
    
    # Rename -- set the name to the cell id
    merged_sigobj._name = cell_id
    
    with open(outsig, 'wt') as f:
        sourmash.save_signatures([merged_sigobj], fp=f)
    return cell_id, True
        

def merge_aligned_unaligned_sigs(sig_folder, outdir_base, verbose=False, n_jobs=32):

    sig_df = get_sig_file_df(sig_folder, verbose=verbose)
    if verbose:
        print("sig_df.head()\n", sig_df.head())

    ## Filter for only signatures that have both aligned and unaligned versions
    # Both aligned and unaligned present
    sig_df_both_aligned_unaligned = sig_df.groupby(SKETCH_CELL_PARAMS).filter(lambda x: len(x) == 2)
    if verbose:
        print('--\nsig_df_both_aligned_unaligned.shape:', sig_df_both_aligned_unaligned.shape)
    sig_df_both_aligned_unaligned.head()

    if verbose:
        print("--\nsig_df_both_aligned_unaligned.cell_id.value_counts():\n", sig_df_both_aligned_unaligned.cell_id.value_counts())

    ## Merge signatures from same channel, cell barcode, ksize, molecule, scaled value
    grouped = sig_df_both_aligned_unaligned.groupby(SKETCH_CELL_PARAMS)

    merged_success = Parallel(n_jobs=n_jobs)(delayed(merge_signatures)(params, df, outdir_base=outdir_base) for params, df in tqdm(grouped))
    merged_success = pd.Series(dict(merged_success))
    return merged_success
