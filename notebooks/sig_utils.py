from collections import Counter
import glob
import itertools
import os
import sys
import logging

from joblib import Parallel, delayed
import pandas as pd
from khmer import Countgraph
import sourmash
from sourmash.logging import error, notify, set_quiet
from sourmash import sourmash_args


from tqdm import tqdm
import sig2kmer

logging.basicConfig(filename="sig_utils.log", level=logging.DEBUG)

KSIZES = 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75
MOLECULES = "protein", "dayhoff"
SCALEDS = 2, 5, 10


SKETCH_INFO_PATTERN = "(?P<mol_or_alpha>molecule|alphabet)-(?P<alphabet>\w+)__ksize-(?P<ksize>\d+)__(?P<sketch_style>num_hashes|scaled)-(?P<sketch_value>\d+)"

SKETCH_PARAMS = ["alphabet", "ksize", "sketch_style", "sketch_value"]

CELL_PARAMS = ["channel", "cell_barcode"]
SKETCH_CELL_PARAMS = SKETCH_PARAMS + CELL_PARAMS


def sanitize(name):
    """Make the name able to be a valid path name.

    No spaces or slashes, everything lowercase"""
    return name.lower().replace(" ", "_").replace("/", "-slash-")


def make_sketch_id(alpha, ksize, style, value):
    sketch_id = f"alphabet-{alpha}__ksize-{ksize}__{style}-{value}"
    return sketch_id


def get_sig_file_df(sig_folder, verbose=False, aligned_unaligned_merged=False):
    df = pd.Series(glob.glob(f"{sig_folder}/*.sig"), name="fullpath").to_frame()
    df["basename"] = df["fullpath"].map(os.path.basename)
    if verbose:
        print("df.shape:", df.shape)

    if aligned_unaligned_merged:
        cell_info = pd.Series(
            df["basename"].str.split(".").str[0], name="cell_id"
        ).to_frame()
        sig_info = pd.DataFrame(index=df.index)
    else:
        sig_info = df["basename"].str.extractall(SKETCH_INFO_PATTERN)
        sig_info = sig_info.droplevel(-1)
        sig_info["ksize"] = sig_info["ksize"].astype(int)
        sig_info["sketch_value"] = sig_info["sketch_value"].astype(int)
        cell_info = df["basename"].str.split("__", expand=True)
        cols_to_drop = cell_info.columns.intersection([3, 4, 5, 6, 7])
        cell_info = cell_info.drop(cols_to_drop, axis=1)
        cell_info = cell_info.rename(
            columns={0: "channel", 1: "alignment_status", 2: "cell_barcode"}
        )
        cell_info["cell_id"] = cell_info.apply(
            lambda x: "{channel}__{cell_barcode}".format(**x), axis=1
        )

    sigs_with_metadata = pd.concat([df, sig_info, cell_info], axis=1)

    # don't worrry about sorting for now
    #     sort_cols = sigs_with_metadata.columns.intersection(SKETCH_CELL_PARAMS)
    #     if len(sort_cols) > 0:
    #         sigs_with_metadata = sigs_with_metadata.sort_values(sort_cols)
    if verbose:
        print("sigs_with_metadata.shape:", sigs_with_metadata.shape)
    return sigs_with_metadata


## -- Sourmash signature stuff --- ###


# from sourmash.sig import _check_abundance_compatibility


def unset_abundances_and_update_bloomfilter(minhash, bloomfilter):
    #     # Set abundances to 1 to only count that the k-mer is present
    #     try:
    #         # Make sure not to edit the underlying object
    #         mh_copy = minhash.copy()
    #         mh_copy.set_abundances({hashval: 1 for hashval in mh_copy.hashes.items()})
    #     except AttributeError:
    #         # Abundances aren't set, so use the original object
    #         mh_copy = minhash

    bloomfilter.update(minhash)


def merge(
    filenames,
    ksize,
    moltype,
    name=None,
    flatten=False,
    outsig=None,
    min_sig_fraction=None,
):
    """
    merge one or more signatures.

    Adapted from 'sourmash sig merge' command


    min_sig_fraction : float or int
        Proportion of k-mer signatures that each hash must be present in to be counted
        If float, then assumed to be a proportion
        If int, then assumed to be the number of cells
    """

    first_sig = None
    mh = None
    total_loaded = 0

    filter_hashes = min_sig_fraction is not None

    if filter_hashes:
        bloomfilter = Countgraph(ksize, starting_size=1e7, n_tables=4)

    for sigfile in filenames:
        this_n = 0

        for sigobj in sourmash.load_file_as_signatures(
            sigfile, ksize=ksize, select_moltype=moltype
        ):

            if sigobj is None:
                error(
                    "No signature in file {}",
                    sigfile,
                )

            # first signature? initialize a bunch of stuff
            if first_sig is None:
                first_sig = sigobj
                mh = first_sig.minhash.copy_and_clear()
            #                 counter = Counter()

            try:
                sigobj_mh = sigobj.minhash

                if filter_hashes:
                    for hashval in sigobj_mh.hashes:
                        bloomfilter.count(hashval)
                #                     unset_abundances_and_update_bloomfilter(sigobj_mh, bloomfilter)
                #                     counter.update(sigobj_mh.hashes.keys())

                mh.merge(sigobj_mh)
            except:
                error(
                    "ERROR when merging signature '{}' ({}) from file {}",
                    sigobj.name(),
                    sigobj.md5sum()[:8],
                    sigfile,
                )
                raise

            this_n += 1
            total_loaded += 1

    merged_sigobj = sourmash.SourmashSignature(mh)

    if filter_hashes:
        if isinstance(min_sig_fraction, float):
            min_observations_per_kmer = min_sig_fraction * total_loaded
        else:
            min_observations_per_kmer = min_sig_fraction
        # Retain only k-mers that have been observed in enough samples
        hashes_to_keep = {
            hashval: abundance
            for hashval, abundance in merged_sigobj.minhash.hashes.items()
            if bloomfilter.get(hashval) >= min_observations_per_kmer
        }
        filtered_mh = merged_sigobj.minhash.copy_and_clear()
        filtered_mh.set_abundances(hashes_to_keep)
        merged_sigobj = sourmash.SourmashSignature(filtered_mh)
    #         merged_sigobj.minhash.remove_many(hashes_to_remove)
    #     import pdb

    #     pdb.set_trace()

    if name is not None:
        merged_sigobj._name = name

    if outsig is not None:
        with open(outsig, "wt") as f:
            sourmash.save_signatures([merged_sigobj], fp=f)

    return merged_sigobj


def merge_signatures(params, df, outdir_base):
    alpha, ksize, style, value, channel, barcode = params

    sketch_id = f"alphabet-{alpha}__ksize-{ksize}__{style}-{value}"
    outdir = f"{outdir_base}/{sketch_id}"
    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except FileExistsError:
            # Already created by another thread
            pass

    cell_id = f"{channel}__{barcode}"
    outsig = f"{outdir}/{cell_id}.sig"
    insigs = " ".join(df["fullpath"])

    try:
        merged_sigobj = merge(df["fullpath"].values, ksize=ksize, moltype=alpha)
    except ValueError:
        print(f"cell: {cell_id} was not able to be merged\n")
        return cell_id, False

    # Rename -- set the name to the cell id
    merged_sigobj._name = cell_id

    with open(outsig, "wt") as f:
        sourmash.save_signatures([merged_sigobj], fp=f)
    return cell_id, True


def merge_aligned_unaligned_sigs(sig_folder, outdir_base, verbose=False, n_jobs=32):

    sig_df = get_sig_file_df(sig_folder, verbose=verbose)
    if verbose:
        print("sig_df.head()\n", sig_df.head())

    ## Filter for only signatures that have both aligned and unaligned versions
    # Both aligned and unaligned present
    sig_df_both_aligned_unaligned = sig_df.groupby(SKETCH_CELL_PARAMS).filter(
        lambda x: len(x) == 2
    )
    if verbose:
        print(
            "--\nsig_df_both_aligned_unaligned.shape:",
            sig_df_both_aligned_unaligned.shape,
        )
    sig_df_both_aligned_unaligned.head()

    if verbose:
        print(
            "--\nsig_df_both_aligned_unaligned.cell_id.value_counts():\n",
            sig_df_both_aligned_unaligned.cell_id.value_counts(),
        )

    ## Merge signatures from same channel, cell barcode, ksize, molecule, scaled value
    grouped = sig_df_both_aligned_unaligned.groupby(SKETCH_CELL_PARAMS)

    merged_success = Parallel(n_jobs=n_jobs)(
        delayed(merge_signatures)(params, df, outdir_base=outdir_base)
        for params, df in tqdm(grouped)
    )
    merged_success = pd.Series(dict(merged_success))
    return merged_success


def _merge_signatures(
    params,
    df,
    outdir_base,
    ksize=None,
    moltype=None,
    style="scaled",
    value=10,
    sig_fullpath="sig_fullpath",
    force=False,
):
    channel, barcode = params

    sketch_id = f"alphabet-{moltype}__ksize-{ksize}__{style}-{value}"
    outdir = f"{outdir_base}/{sketch_id}"
    cell_id = f"{channel}__{barcode}"
    outsig = f"{outdir}/{cell_id}.sig"

    if not force and os.path.exists(outsig):
        return cell_id, True

    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except FileExistsError:
            # Already created by another thread
            pass

    #     insigs = ' '.join(df['sig_fullpath'])

    try:
        merged_sigobj = merge(df[sig_fullpath].values, ksize=ksize, moltype=moltype)
    except ValueError:
        print(f"cell: {cell_id} was not able to be merged\n")
        return cell_id, False

    # Rename -- set the name to the cell id
    merged_sigobj._name = cell_id

    with open(outsig, "wt") as f:
        sourmash.save_signatures([merged_sigobj], fp=f)
    return cell_id, True


def _flatten(
    filename, ksize, moltype, cell_id, md5=None, name=None, quiet=True  # to rename
):
    set_quiet(quiet)

    siglist = [filename]
    outlist = []
    total_loaded = 0
    siglist = list(siglist)
    # raise ValueError
    total_loaded += len(siglist)

    # select!
    if md5 is not None:
        siglist = [ss for ss in siglist if md5 in ss.md5sum]

    if name is not None:
        siglist = [ss for ss in siglist if name in ss.name]

    for ss in siglist:
        flattened_mh = ss.minhash.copy_and_clear()
        flattened_mh.track_abundance = False
        flattened_mh.add_many(ss.minhash.get_mins())

        ss.minhash = flattened_mh

    outlist.extend(siglist)
    flattened = sourmash.save_signatures(outlist)

    flattened = sourmash.load_one_signature(
        flattened, ksize=ksize, select_moltype=moltype
    )
    flattened._name = cell_id

    return flattened


def _subtract(
    original_sig, sig_to_subtract_path, ksize, cell_id, moltype="dayhoff", quiet=True
):
    """
    subtract one or more signatures from another
    """
    set_quiet(quiet)
    from_sigobj = original_sig

    from_mh = from_sigobj.minhash
    subtract_mins = set(from_mh.hashes)

    total_loaded = 0

    if isinstance(sig_to_subtract_path, sourmash.SourmashSignature):
        sigobj = sig_to_subtract_path
        subtract_mins -= set(sigobj.minhash.hashes)
        total_loaded += 1
        logging.info("loaded and subtracted signatures from {}...".format(sigobj.name))
    elif isinstance(sig_to_subtract_path, set):
        subtract_mins -= sig_to_subtract_path
        total_loaded += 1
    else:
        for sigobj in sourmash_args.load_file_as_signatures(
            sig_to_subtract_path, ksize=ksize, select_moltype=moltype
        ):
            if not sigobj.minhash.is_compatible(from_mh):
                logging.debug(
                    f"incompatible minhashes from {sig_to_subtract_path}; specify -k and/or molecule type."
                )
                return

            total_loaded += 1
        logging.info("loaded and subtracted signatures from {}...".format(sigobj.name))

    if not total_loaded:
        logging.debug(
            f"no signatures to subtract from {original_sig} {sig_to_subtract_path}t!?"
        )
        return

    subtract_mh = from_sigobj.minhash.copy_and_clear()
    subtract_mh.add_many(subtract_mins)

    subtract_sigobj = sourmash.SourmashSignature(subtract_mh)
    subtract_sigobj._name = cell_id

    return subtract_sigobj


def _intersect(
    ksize,
    signatures_to_merge: list,
    abund_sig,
    cell_id,  # to rename
    moltype="dayhoff",
    quiet=True,
):
    """
    intersect one or more signatures by taking the intersection of hashes.
    This function always removes abundances.
    """
    set_quiet(quiet)

    first_sig = None
    mins = None
    total_loaded = 0

    for sigobj in signatures_to_merge:
        if first_sig is None:
            first_sig = sigobj
            mins = set(sigobj.minhash.hashes)
        else:

            # check signature compatibility --
            if not sigobj.minhash.is_compatible(first_sig.minhash):
                raise ValueError(
                    "incompatible minhashes; specify -k and/or molecule type."
                )

        mins.intersection_update(sigobj.minhash.hashes)
        total_loaded += 1
        logging.debug(f"loaded and intersected signatures from {sigobj.name}...")

    if total_loaded == 0:
        raise ValueError("no signatures to merge!?")

    if not abund_sig.minhash.track_abundance:
        logging.debug(f" abundances not set for: ksize {ksize} and cell_id {cell_id}")
        pass
        # raise ValueError("--track-abundance not set on loaded signature?! exiting.")
    else:
        intersect_mh = abund_sig.minhash.copy_and_clear()
        abund_mins = abund_sig.minhash.hashes

        # do one last intersection
        mins.intersection_update(abund_mins)
        abund_mins = {k: abund_mins[k] for k in mins}

        intersect_mh.set_abundances(abund_mins)
        intersect_sigobj = sourmash.SourmashSignature(intersect_mh)
        intersect_sigobj._name = cell_id
        logging.debug(f" intersect passed for: ksize {ksize}, cell_id {cell_id}")
        return intersect_sigobj


def parallel_merge_aligned_unaligned_sigs(
    sig_df,
    outdir_base,
    groupby=["channel", "cell_barcode"],
    verbose=False,
    n_jobs=32,
    ksizes=KSIZES,
    moltype="dayhoff",
    sig_fullpath="sig_fullpath",
    force=False,
):

    ## Merge signatures from same channel, cell barcode, ksize, molecule, scaled value
    grouped = sig_df.groupby(groupby)
    n_iter = len(grouped) * len(ksizes)
    print(f"total number of iterations: {n_iter}")

    merged_success = Parallel(n_jobs=n_jobs)(
        delayed(_merge_signatures)(
            params,
            df,
            outdir_base=outdir_base,
            ksize=ksize,
            moltype=moltype,
            sig_fullpath=sig_fullpath,
        )
        for (params, df), ksize in tqdm(
            itertools.product(grouped, ksizes), total=n_iter
        )
    )
    merged_success = pd.Series(dict(merged_success))
    return merged_success


def get_hashes_to_remove(sigs, percent_threshold, ksize):
    shared_hashes = Countgraph(ksize, starting_size=1e7, n_tables=4)
    n_sigs_threshold = percent_threshold * len(sigs)
    hashes_to_remove = set([])

    for sig in sigs:
        for hashval in sig.minhash.hashes:
            shared_hashes.count(hashval)
            if shared_hashes.get(hashval) >= n_sigs_threshold:
                hashes_to_remove.add(hashval)

    return hashes_to_remove


def remove_hashes(sig, hashes_to_remove):
    subtract_mh = sig.minhash.copy_and_clear()
    hash_dict = sig.minhash.hashes
    hashes_to_keep = set(hash_dict.keys()) - hashes_to_remove
    subtract_mh.add_many(hashes_to_keep)
    subtract_mh.set_abundances({h: hash_dict[h] for h in hashes_to_keep})
    subtract_sigobj = sourmash.SourmashSignature(subtract_mh, name=sig.name())
    return subtract_sigobj


def remove_common_hashes(sigs, percent_threshold, ksize):
    hashes_to_remove = get_hashes_to_remove(sigs, percent_threshold, ksize)
    subtracted_sigs = [remove_hashes(sig, hashes_to_remove) for sig in sigs]

    return subtracted_sigs


def load_sigfiles(sigfiles, ksize=None, moltype=None):
    sigs = []
    for filename in sigfiles:
        sigs.extend(
            sourmash.load_file_as_signatures(
                filename, ksize=ksize, select_moltype=moltype
            )
        )
    return sigs


def remove_common_hashes_from_sig_df(
    sig_df,
    sketch_id,
    ksize,
    moltype,
    fraction_threshold=0.8,
    sig_col="sig_path",
    output_dir=None,
    force=False,
    output_subfolder=None,
    #     create_hash_count_csv=True,
    use_sparse_matrix=True,
):
    """
    Remove hashes that are present in `fraction_threshold` of all signatures, e.g in 80% or more of signatures

    sig_df : pandas.DataFrame
        A pandas dataframe containing the filename paths for merged cell type signatures
    sketch_id : str
        A string of format f"alphabet-{alpha}__ksize-{ksize}__{style}-{value}"
    ksize : int
        The k-mer size of interest
    output_dir : str, optional
        If set, then write the signatures to:
        {output_dir}/{sketch_id}/{sig_df[sig_col].map.os.path(basename)}.sig
    output_subfolder : str, optional
        IF set and output_dir is set, then write the signatures to:
        {output_dir}/{output_subfolder}/{sketch_id}/{sig_df[sig_col].map.os.path(basename)}.sig
    """
    sigs = load_sigfiles(sig_df[sig_col].values, ksize, moltype)
    sigs_without_common_hashes = remove_common_hashes(sigs, fraction_threshold)
    series = pd.Series(
        sigs_without_common_hashes,
        index=[sig.name() for sig in sigs]
    )
    if output_dir:
        if output_subfolder is None:
            sketch_dir = os.path.join(output_dir, sketch_id)
        else:
            sketch_dir = os.path.join(output_dir, output_subfolder, sketch_id)
        print(f"Writing to {sketch_dir}")
        # Early exit if force=False
        if os.path.exists(sketch_dir) and not force:
            return

        if not os.path.exists(sketch_dir):
            try:
                os.makedirs(sketch_dir)
            except FileExistsError:
                # Some other thread made this
                pass

    sigs = load_sigfiles(sig_df[sig_col].values, ksize, moltype)

    sigs_without_common_hashes = remove_common_hashes(sigs, fraction_threshold, ksize)
    series = pd.Series(sigs_without_common_hashes, index=[sig.name() for sig in sigs])

    if output_dir:
        basenames = sig_df[sig_col].map(os.path.basename)
        for subtracted_sig, basename in zip(sigs_without_common_hashes, basenames):
            outfile = os.path.join(sketch_dir, basename)
            with open(outfile, "wt") as f:
                sourmash.save_signatures([subtracted_sig], fp=f)

    return sigs_without_common_hashes


def extract_channel_alignment_barcode(singlecell_fasta, double_aligned=False):
    basename = os.path.basename(singlecell_fasta)
    if double_aligned:
        # "__aligned__aligned__" or "__unaligned__unaligned__"
        channel, alignment_status, alignment_status2, barcode, suffix = basename.split(
            "__"
        )
    else:
        # "__aligned__" or "__unaligned__"
        channel, alignment_status, barcode, suffix = basename.split("__")
    return channel.strip("_"), alignment_status.strip("_"), barcode.strip("_")


def extract_cell_id_from_fasta(singlecell_fasta, double_aligned=False):
    channel, alignment_status, barcode = extract_channel_alignment_barcode(
        singlecell_fasta, double_aligned=double_aligned
    )
    return f"{channel}__{barcode}"


def clean_fasta_name(
    basename,
    strings_to_remove=[
        "__aligned",
        "__possorted_genome_bam",
        "_possorted_genome_bam",
        "__unaligned",
    ],
):
    """Remove alignment status and extra bam text from orpheum-created translated coding reads"""
    new_name = None
    for to_remove in strings_to_remove:
        if new_name is None:
            # First time --> take original basename
            new_name = basename.replace(to_remove, "")
        else:
            new_name = new_name.replace(to_remove, "")

    new_name = new_name.split("_coding_reads")[0].strip("_")
    return new_name


def get_hashvals_in_singlecell_fasta(
    singlecell_fasta,
    cell_ontology_classes,
    query_hashvals,
    seqout_dir,
    ksize,
    moltype,
    gene_name_tag="GN",
    double_aligned=False,
):
    channel, alignment_status, barcode = extract_channel_alignment_barcode(
        singlecell_fasta, double_aligned=double_aligned
    )

    cell_id = f"{channel}__{barcode}"

    # Only consider cells with annotation
    if cell_id not in cell_ontology_classes:
        logging.info(f"cell id {cell_id} is not annotated")
        return

    cell_ontology = sanitize(cell_ontology_classes[cell_id])

    seqout_filename = os.path.join(
        seqout_dir,
        alignment_status,
        cell_ontology,
        f"{cell_ontology}---{cell_id}.fasta",
    )
    dirname = os.path.dirname(seqout_filename)

    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(seqout_filename, "wt") as f:
        df = sig2kmer.get_kmers_in_seqfiles(
            [singlecell_fasta],
            query_hashvals,
            ksize=ksize,
            moltype=moltype,
            seqout_fp=f,
        )
    df["gene_name"] = df["read_name"].str.extract(f"{gene_name_tag}:Z:([\w-]+)")
    df["cell_id"] = cell_id
    df["alignment_status"] = alignment_status
    return df


def merge_celltype_sigs(
    sketch_id,
    ksize,
    moltype,
    celltype_name,
    df,
    merged_celltype_outdir_base,
    force=False,
    dryrun=False,
    sig_path_col="sig_path",
    min_sig_fraction=None,
):
    celltype_name_sanitized = sanitize(celltype_name)
    sketch_id_outdir = os.path.join(merged_celltype_outdir_base, sketch_id)

    output_sig = os.path.join(sketch_id_outdir, celltype_name_sanitized + ".sig")
    if dryrun:
        print(f"Merging cell signatures to:\n{output_sig}")
        return output_sig

    if not os.path.exists(sketch_id_outdir):
        try:
            os.makedirs(sketch_id_outdir)
        except FileExistsError:
            # Some other thread made this in between detecting the
            # folder didn't exist and actually making it
            pass
    if force or (not os.path.exists(output_sig) and not force):
        filenames = df[sig_path_col]
        #         import pdb; pdb.set_trace()
        merge(
            filenames,
            ksize=ksize,
            moltype=moltype,
            outsig=output_sig,
            name=celltype_name,
            min_sig_fraction=min_sig_fraction,
        )
    return output_sig
