import csv
import itertools
import os
import sys

from joblib import Parallel, delayed
import pandas as pd
from sourmash import MinHash, load_sbt_index, create_sbt_index
from sourmash import signature as sig
from sourmash import sourmash_args
from sourmash.logging import notify, error, print_results, set_quiet
from sourmash.search import search_databases
from sourmash.index import LinearIndex
from tqdm import tqdm

import sig_utils

KSIZES = 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75
MOLECULES = "protein", "dayhoff"
SCALEDS = 2, 5, 10

CELLTYPE_COLS = [
    "cell_ontology_class",
    "narrow_group",
    "broad_group",
    "compartment_group",
]

MATCHING_CELLTYPE_COLS = [f"{x}_matches" for x in CELLTYPE_COLS]
CELLTYPE_TO_MATCHING_COLS = dict(zip(CELLTYPE_COLS, MATCHING_CELLTYPE_COLS))


SHARED_CELLTYPES = set(['Alveolar Epithelial Type 2',
 'B cell',
 'Capillary',
 'Dendritic',
 'Fibroblast',
 'Macrophage',
 'Monocyte',
 'Natural Killer T cell',
 'Smooth Muscle and Myofibroblast',
 'T cell'])
SHARED_CELLTYPES

def calculate_moltype(dna=False, dayhoff=False, hp=False, protein=False, default=None):
    moltype = default

    n = 0
    if dna:
        moltype = "DNA"
        n += 1
    if dayhoff:
        moltype = "dayhoff"
        n += 1
    if hp:
        moltype = "hp"
        n += 1
    if protein:
        moltype = "protein"
        n += 1

    if n > 1:
        error("cannot specify more than one of --dna/--rna/--protein/--hp/--dayhoff")
        sys.exit(-1)

    return moltype


def load_sigs_as_index(filenames, moltype, ksize):
    n_signatures = 0
    n_databases = 0
    databases = []

    for filename in filenames:
        notify("loading from {}...", filename, end="\r")

        try:
            db, dbtype = sourmash_args._load_database(
                filename, traverse=False, traverse_yield_all=False
            )
        except IOError as e:
            notify(str(e))
            sys.exit(-1)

        # are we collecting signatures from a directory/path?
        # NOTE: error messages about loading will now be attributed to
        # directory, not individual file.
        if os.path.isdir(filename):
            assert dbtype == sourmash_args.DatabaseType.SIGLIST

            siglist = sourmash_args._select_sigs(db, moltype=moltype, ksize=ksize)
            # siglist = sourmash_args.filter_compatible_signatures(query, siglist, 1)
            linear = LinearIndex(siglist, filename=filename)
            databases.append((linear, filename, False))

            n_signatures += len(linear)

        if dbtype == sourmash_args.DatabaseType.SIGLIST:
            siglist = sourmash_args._select_sigs(db, moltype=moltype, ksize=ksize)
            # siglist = sourmash_args.filter_compatible_signatures(query, siglist, False)
            siglist = list(siglist)

            if not siglist:
                notify("--\nno compatible signatures found in:\n'{}'", filename)
                for sig in siglist:
                    notify(
                        "molecule: {}, ksize: {}",
                        sig.minhash.moltype,
                        sig.minhash.ksize,
                    )
                notify("--")
                # Don't exit if only one file has incompatible signatures
            #                 sys.exit(-1)

            linear = LinearIndex(siglist, filename=filename)
            databases.append((linear, filename, "signature"))

            notify("loaded {} signatures from {}", len(linear), filename, end="\r")
            n_signatures += len(linear)

    notify(" " * 79, end="\r")
    if n_signatures and n_databases:
        notify(
            "loaded {} signatures and {} databases total.", n_signatures, n_databases
        )
    elif n_signatures:
        notify("loaded {} signatures.", n_signatures)
    elif n_databases:
        notify("loaded {} databases.", n_databases)
    else:
        notify("** ERROR: no signatures or databases loaded?")
        sys.exit(-1)

    if databases:
        print("")

    return databases


def search_singletons(
    queryfile,
    databases,
    ksize,
    dna=False,
    dayhoff=False,
    hp=False,
    protein=False,
    md5=None,
    scaled=0,
    containment=False,
    save_matches=None,
    output=None,
    ignore_abundance=False,
    threshold=0.08,
    best_only=False,
    num_results=3,
    quiet=False,
    use_tqdm=False,
):
    set_quiet(quiet)
    moltype = calculate_moltype(dna, dayhoff, hp, protein)
    notify("Moltype: {}", moltype)
    notify("Ksize: {}", ksize)

    # set up the queries
    queries = sourmash_args.load_file_as_signatures(
        queryfile,
        ksize=ksize,
        select_moltype=moltype,  # select_md5=md5
    )
    notify(
        "loaded queries",
    )

    # set up the search databases
    databases = load_sigs_as_index(databases, moltype=moltype, ksize=ksize)
    notify("Number of database signatures loaded: {}", len(databases))

    if not len(databases):
        error("Nothing found to search!")
        sys.exit(-1)

    # do the actual search
    results = {}

    if use_tqdm:
        query_iterator = tqdm(queries)
    else:
        query_iterator = queries
    for query in query_iterator:
        print_results("\n- Query: {} -", query.name())
        results[query.name()] = search_query(
            query,
            databases,
            threshold,
            containment,
            best_only,
            ignore_abundance,
            output,
            save_matches,
            num_results,
            scaled=scaled,
        )

    results_df = pd.DataFrame(
        {
            key: {sr.match.name(): sr.similarity * 100 for sr in search_results}
            for key, search_results in results.items()
        }
    )
    return results_df


def search_query(
    query,
    databases,
    threshold,
    containment,
    best_only,
    ignore_abundance,
    output,
    save_matches,
    num_results,
    scaled=0,
):
    # downsample if requested
    if scaled:
        if query.minhash.max_hash == 0:
            error("cannot downsample a signature not created with --scaled")
            sys.exit(-1)

        if scaled != query.minhash.scaled:
            notify(
                "downsampling query from scaled={} to {}",
                query.minhash.scaled,
                int(scaled),
            )
        query.minhash = query.minhash.downsample(scaled=scaled)

    # forcibly ignore abundances if query has no abundances
    if not query.minhash.track_abundance:
        ignore_abundance = True

    # assume compatibility for now..
    # sourmash_args.check_tree_is_compatible("treename", databases, query,
    #                                        is_similarity_query=not containment)

    #     import pdb; pdb.set_trace()
    results = search_databases(
        query,
        databases,
        threshold,
        containment,
        best_only,
        ignore_abundance,
        unload_data=True,
    )

    n_matches = len(results)
    if best_only:
        num_results = 1

    if not num_results or n_matches <= num_results:
        print_results("{} matches:".format(len(results)))
    else:
        print_results("{} matches; showing first {}:", len(results), num_results)
        n_matches = num_results

    # output!
    print_results("similarity   match")
    print_results("----------   -----")
    for sr in results[:n_matches]:
        pct = "{:.1f}%".format(sr.similarity * 100)
        name = sr.match._display_name(60)
        print_results("{:>6}       {}", pct, name)

    if best_only:
        notify("** reporting only one match because --best-only was set")

    if output:
        fieldnames = ["similarity", "name", "filename", "md5"]

        with sourmash_args.FileOutput(output, "wt") as fp:
            w = csv.DictWriter(fp, fieldnames=fieldnames)

            w.writeheader()
            for sr in results:
                d = dict(sr._asdict())
                del d["match"]
                w.writerow(d)

    # save matching signatures upon request
    if save_matches:
        notify('saving all matched signatures to "{}"', save_matches)
        with sourmash_args.FileOutput(save_matches, "wt") as fp:
            sig.save_signatures([sr.match for sr in results], fp)
    return results


def load_matching_signatures_into_tree(
    filenames,
    ksize,
    moltype,
    scaled=0,
    append=False,
    sbt_name=None,
    bf_size=1e5,
    n_children=2,
    return_n=False,
):
    if append:
        tree = load_sbt_index(sbt_name)
    else:
        tree = create_sbt_index(bf_size, n_children=n_children)

    n = 0
    ksizes = set()
    moltypes = set()
    nums = set()
    scaleds = set()
    for f in filenames:
        if n % 100 == 0:
            notify("\r...reading from {} ({} signatures so far)", f, n, end="")
        siglist = sig.load_signatures(f, ksize=ksize, select_moltype=moltype)

        # load all matching signatures in this file
        ss = None
        for ss in siglist:
            ksizes.add(ss.minhash.ksize)
            moltypes.add(sourmash_args.get_moltype(ss))
            nums.add(ss.minhash.num)

            if scaled:
                ss.minhash = ss.minhash.downsample_scaled(scaled)
            scaleds.add(ss.minhash.scaled)

            tree.insert(ss)
            n += 1

        if not ss:
            continue

        check_signature_compatibilty_to_tree(ksizes, moltypes, nums, scaleds)
    notify("")

    # did we load any!?
    if n == 0:
        error("no signatures found to load into tree!? failing.")
        sys.exit(-1)
    if return_n:
        return tree, n
    else:
        return tree


def check_signature_compatibilty_to_tree(ksizes, moltypes, nums, scaleds):
    # check to make sure we aren't loading incompatible signatures
    if len(ksizes) > 1 or len(moltypes) > 1:
        error("multiple k-mer sizes or molecule types present; fail.")
        error("specify --dna/--protein and --ksize as necessary")
        error(
            "ksizes: {}; moltypes: {}", ", ".join(map(str, ksizes)), ", ".join(moltypes)
        )
        sys.exit(-1)
    if nums == {0} and len(scaleds) == 1:
        pass  # good
    elif scaleds == {0} and len(nums) == 1:
        pass  # also good
    else:
        error("trying to build an SBT with incompatible signatures.")
        error("nums = {}; scaleds = {}", repr(nums), repr(scaleds))
        sys.exit(-1)


def do_sourmash_search(query_sig, database_sigs, ksize, molecule, scaled=5):
    protein = False
    dayhoff = False
    if molecule == "protein":
        protein = True
        dayhoff = False
    elif molecule == "dayhoff":
        protein = False
        dayhoff = True
    else:
        raise ValueError(
            f"Unrecognized molecule provided: '{molecule}', must be one of 'dayhoff' or 'protein'"
        )

    df = sourmash_search_utils.search_singletons(
        query_sig,
        database_sigs,
        ksize=ksize,
        dayhoff=dayhoff,
        protein=protein,
        containment=True,
        quiet=True,
        scaled=scaled,
    )
    if not df.empty:
        tidy = df.unstack().reset_index()
        tidy = tidy.rename(
            columns={"level_0": "refseq_query", "level_1": "fasta_id", 0: "containment"}
        )
        fasta_info = tidy["fasta_id"].str.split("__", expand=True)
        tidy["channel"] = fasta_info[0]
        tidy["alignment_status"] = fasta_info[1]
        tidy["cell_barcode"] = fasta_info[2]
        tidy["ksize"] = ksize
        tidy["molecule"] = molecule
        tidy["scaled"] = scaled
        tidy["cell_id"] = tidy.apply(
            lambda x: "{channel}__{cell_barcode}".format(**x), axis=1
        )
        return tidy


def do_sourmash_search(query_sig, database_sigs, ksize, molecule, scaled=5, **kwargs):
    if molecule == "protein":
        protein = True
        dayhoff = False
    elif molecule == "dayhoff":
        protein = False
        dayhoff = True

    df = search_singletons(
        query_sig,
        database_sigs,
        ksize=ksize,
        dayhoff=dayhoff,
        protein=protein,
        containment=True,
        quiet=True,
        scaled=scaled,
        **kwargs,
    )
    if not df.empty:
        tidy = df.unstack().reset_index()
        tidy = tidy.rename(
            columns={"level_0": "refseq_query", "level_1": "fasta_id", 0: "containment"}
        )
        fasta_info = tidy["fasta_id"].str.split("__", expand=True)
        tidy["channel"] = fasta_info[0]
        tidy["alignment_status"] = fasta_info[1]
        tidy["cell_barcode"] = fasta_info[2]
        tidy["ksize"] = ksize
        tidy["molecule"] = molecule
        tidy["scaled"] = scaled
        tidy["cell_id"] = tidy.apply(
            lambda x: "{channel}__{cell_barcode}".format(**x), axis=1
        )
        return tidy


def search_parallel(
    query_sig,
    database_sigs,
    cell_annotation,
    ksizes=KSIZES,
    molecules=MOLECULES,
    scaleds=SCALEDS,
    n_jobs=64,
    **kwargs,
):
    dfs = Parallel(n_jobs=n_jobs)(
        delayed(do_sourmash_search)(
            query_sig, database_sigs, ksize, molecule, scaled, **kwargs
        )
        for ksize, molecule, scaled in tqdm(
            itertools.product(ksizes, molecules, scaleds)
        )
    )

    results = pd.concat(dfs)

    # add cell type annotation
    results_annotated = results.join(cell_annotation, rsuffix="_x", on="cell_id")

    return results_annotated


def add_is_matching_col(
    df,
    renamer1,
    renamer2,
    celltype_to_matching_cols=CELLTYPE_TO_MATCHING_COLS,
):
    ## Add matching col of whether the category matches
    for celltype_col, matching_col in celltype_to_matching_cols.items():
        celltype1_col = renamer1[celltype_col]
        celltype2_col = renamer2[celltype_col]
        df[matching_col] = df[celltype1_col] == df[celltype2_col]
    return df


def read_and_process_sourmash_search_results(
    csv,
    query_species,
    adata,
    column_renamer,
    shared_celltypes=SHARED_CELLTYPES,
    celltype_col="broad_group",
    sbt_species="mouse",
    celltype_cols=CELLTYPE_COLS,
):
    query_cell_id_col = f"{query_species}_cell_id"
    query_celltype_col = f"{query_species}_{celltype_col}"

    sbt_celltype_col = f"{sbt_species}_{celltype_col}"

    sketch_id = os.path.basename(os.path.dirname(csv))
    cell_id = os.path.basename(csv).split(".csv")[0]

    df = pd.read_csv(csv)
    if 'containment' in csv:
        df["search_method"] = "containment"
    else:
        df["search_method"] = "similarity_abundance"

    df["sketch_id"] = sketch_id
    sig_info = df["sketch_id"].str.extractall(sig_utils.SKETCH_INFO_PATTERN)
    sig_info = sig_info.droplevel(-1)
    sig_info["ksize"] = sig_info["ksize"].astype(int)
    sig_info["sketch_value"] = sig_info["sketch_value"].astype(int)
    df[query_cell_id_col] = cell_id
    df_metadata = pd.concat([df, sig_info], axis=1)
    df_metadata = df_metadata.drop(["filename", "md5"], axis=1)
    df_metadata = df_metadata.rename(columns={"name": sbt_celltype_col})
    df_metadata[f"{query_species}_channel"] = (
        df_metadata[query_cell_id_col].str.split("__").str[0]
    )

    # Add cell_ontology_class, narrow_group, broad_group, compartment_group for each mouse
    df_metadata_ref_query_ontologies = add_ontologies(
        df_metadata, adata, query_cell_id_col, column_renamer, celltype_cols
    )
    if shared_celltypes is not None:
        df_metadata_ref_query_ontologies = df_metadata_ref_query_ontologies.query(
            f"{sbt_celltype_col} in @SHARED_CELLTYPES and {query_celltype_col} in @SHARED_CELLTYPES"
        )
    return df_metadata_ref_query_ontologies


def add_ontologies(
    search_results_df, adata, mouse_cell_id_col, renamer, celltype_cols=CELLTYPE_COLS
):
    search_results_with_ontologies = search_results_df.join(
        adata.obs[celltype_cols], on=mouse_cell_id_col, how="inner"
    )
    # Remove null ksizes
    search_results_with_ontologies["ksize"] = search_results_with_ontologies[
        "ksize"
    ].astype(int)
    search_results_with_ontologies["sketch_value"] = search_results_with_ontologies[
        "sketch_value"
    ].astype(int)
    search_results_with_ontologies.rename(columns=renamer, inplace=True)
    search_results_with_ontologies.reset_index(inplace=True, drop=True)

    return search_results_with_ontologies
