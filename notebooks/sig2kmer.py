#! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

Cribbed from https://github.com/dib-lab/sourmash/pull/724/
"""
import logging
import os
import sys
import argparse
import sourmash

from sourmash.minhash import hash_murmur
import screed
import csv
from sourmash.logging import notify, error, set_quiet
from sourmash.cli.utils import add_construct_moltype_args, add_ksize_arg
from sourmash.sourmash_args import calculate_moltype

from orpheum.sequence_encodings import encode_peptide, AMINO_ACID_SINGLE_LETTERS

import sig_utils
import pandas as pd

NOTIFY_EVERY_BP = 1e7

# logging.basicConfig(filename='sig2kmer.log', level=logging.INFO)


def get_kmer_moltype(sequence, start, ksize, moltype, input_is_protein):
    kmer_in_seq = sequence[start : start + ksize]
    if moltype == "DNA":
        # Get reverse complement
        kmer_rc = screed.rc(kmer_in_seq)
        if kmer_in_seq > kmer_rc:  # choose fwd or rc
            kmer_encoded = kmer_rc
        else:
            kmer_encoded = kmer_in_seq
    elif input_is_protein:
        kmer_encoded = encode_peptide(kmer_in_seq, moltype)
    elif not input_is_protein:
        raise NotImplementedError("Currently cannot translate DNA to protein sequence")
    return kmer_encoded, kmer_in_seq


def revise_ksize(ksize, moltype, input_is_protein):
    """If input is protein, then divide the ksize by three"""
    if moltype == "DNA":
        return ksize
    elif input_is_protein:
        # Ksize includes codons
        return int(ksize / 3)
    else:
        return ksize


def get_kmers_for_hashvals(sequence, hashvals, ksize, moltype, input_is_protein):
    "Return k-mers from 'sequence' that yield hashes in 'hashvals'."
    # uppercase!
    sequence = sequence.upper()

    # Divide ksize by 3 if sequence is protein
    ksize = revise_ksize(ksize, moltype, input_is_protein)

    for start in range(0, len(sequence) - ksize + 1):
        # Skip protein sequences with invalid input
        # (workaround for sencha bug that wrote "Writing translate
        # summary to coding_summary.json" to standard output and thus to the
        # protein fasta)
        if input_is_protein:
            if not all(x in AMINO_ACID_SINGLE_LETTERS for x in sequence):
                continue

        kmer_encoded, kmer_in_seq = get_kmer_moltype(
            sequence, start, ksize, moltype, input_is_protein
        )

        # NOTE: we do not avoid non-ACGT characters, because those k-mers,
        # when hashed, shouldn't match anything that sourmash outputs.
        hashval = hash_murmur(kmer_encoded)
        if hashval in hashvals:
            yield kmer_encoded, kmer_in_seq, hashval


def get_matching_hashes_in_file(
    filename,
    ksize,
    moltype,
    input_is_protein,
    hashes,
    found_kmers,
    m,
    n,
    n_seq,
    seqout_fp,
    watermark,
    first=False,
):
    for record in screed.open(filename):
        n += len(record.sequence)
        n_seq += 1
        while n >= watermark:
            sys.stderr.write("... {} {} {}\r".format(n_seq, watermark, filename))
            watermark += NOTIFY_EVERY_BP

        # now do the hard work of finding the matching k-mers!
        for kmer_encoded, kmer_in_seq, hashval in get_kmers_for_hashvals(
            record.sequence, hashes, ksize, moltype, input_is_protein
        ):
            found_kmers.append([kmer_in_seq, kmer_encoded, hashval, record["name"]])

            # write out sequence
            if seqout_fp:
                seqout_fp.write(
                    ">{}|hashval:{}|kmer:{}|kmer_encoded:{}\n{}\n".format(
                        record.name, hashval, kmer_in_seq, kmer_encoded, record.sequence
                    )
                )
                m += len(record.sequence)
            if first:
                return m, n
    return m, n


def main():
    p = argparse.ArgumentParser()
    p.add_argument("query")  # signature file
    p.add_argument(
        "seqfiles", nargs="+"
    )  # sequence files from which to look for matches
    p.add_argument(
        "--output-sequences",
        type=str,
        default=None,
        help="save matching sequences to this file.",
    )
    p.add_argument(
        "--output-kmers",
        type=str,
        default=None,
        help="save matching kmers to this file.",
    )
    p.add_argument(
        "--input-is-protein",
        action="store_true",
        help="Consume protein sequences - no translation needed.",
    )
    p.add_argument(
        '--quiet',
        action="store_true",
        help='suppress non-error output'
    )
    add_ksize_arg(p)
    add_construct_moltype_args(p)
    args = p.parse_args()

    set_quiet(args.quiet)

    # set up the outputs.
    seqout_fp = None
    if args.output_sequences:
        seqout_fp = open(args.output_sequences, "wt")

    kmerout_fp = None
    if args.output_kmers:
        kmerout_fp = open(args.output_kmers, "wt")
        kmerout_w = csv.writer(kmerout_fp)
        kmerout_w.writerow(
            ["kmer_in_sequence", "kmer_in_alphabet", "hashval", "read_name"]
        )

    # Ensure that protein ksizes are divisible by 3
    if (args.protein or args.dayhoff or args.hp) and not args.input_is_protein:
        if args.ksize % 3 != 0:
            error("protein ksizes must be divisible by 3, sorry!")
            error("bad ksizes: {}", ", ".join(args.ksize))
            sys.exit(-1)

    if not (seqout_fp or kmerout_fp):
        error("No output options given!")
        return -1

    # first, load the signature and extract the hashvals
    moltype = calculate_moltype(args)
    sigobj = sourmash.load_one_signature(
        args.query, ksize=args.ksize, select_moltype=moltype
    )
    query_hashvals = set(sigobj.minhash.hashes)
    query_ksize = sigobj.minhash.ksize

    # now, iterate over the input sequences and output those that overlap
    # with hashes!
    n_seq = 0
    n = 0  # bp loaded
    m = 0  # bp in found sequences
    p = 0  # number of k-mers found
    found_kmers = []
    watermark = NOTIFY_EVERY_BP
    for filename in args.seqfiles:
        m, n = get_matching_hashes_in_file(
            filename,
            query_ksize,
            moltype,
            args.input_is_protein,
            query_hashvals,
            found_kmers,
            m,
            n,
            n_seq,
            seqout_fp,
            watermark,
        )

    if seqout_fp:
        notify("read {} bp, wrote {} bp in matching sequences", n, m)

    if kmerout_fp and found_kmers:
        for kmer_in_seq, kmer_encoded, hashval, read_id in found_kmers:
            kmerout_w.writerow([kmer_in_seq, kmer_encoded, str(hashval), read_id])
        notify("read {} bp, found {} kmers matching hashvals", n, len(found_kmers))


def overlap(sig1, sig2, ksize, moltype):
    """
    provide detailed comparison of two signatures

    from:
    https://github.com/dib-lab/sourmash/blob/ca201cfc1824900bfb94aaadd94e163c70b7cf4e/src/sourmash/sig/__main__.py#L255
    """

    similarity = sig1.similarity(sig2)
    similarity_ignore_abundance = sig1.similarity(sig2, ignore_abundance=True)

    cont1 = sig1.contained_by(sig2)
    cont2 = sig2.contained_by(sig1)

    name1 = sig1.name()
    name2 = sig2.name()

    md5_1 = sig1.md5sum()
    md5_2 = sig2.md5sum()

    ksize = sig1.minhash.ksize
    moltype = sig1.minhash.moltype

    num = sig1.minhash.num
    size1 = len(sig1.minhash)
    size2 = len(sig2.minhash)

    scaled = sig1.minhash.scaled

    hashes_1 = set(sig1.minhash.hashes)
    hashes_2 = set(sig2.minhash.hashes)

    num_common = len(hashes_1.intersection(hashes_2))
    disjoint_1 = len(hashes_1 - hashes_2)
    disjoint_2 = len(hashes_2 - hashes_1)
    num_union = len(hashes_1.union(hashes_2))

    attributes = dict(
        similarity=similarity,
        similarity_ignore_abundance=similarity_ignore_abundance,
        first_contained_in_second=cont1,
        second_contained_in_first=cont2,
        first_n_hashes=size1,
        second_n_hashes=size2,
        n_hashes_in_common=num_common,
        n_hashes_only_in_first=disjoint_1,
        n_hashes_only_in_second=disjoint_2,
        n_hashes_total=num_union,
        first=name1,
        second=name2,
        moltype=moltype,
        ksize=ksize,
    )
    return pd.Series(attributes)


def get_intersecting_hashes(sig1, sig2, abundances_from=None):
    """
    Adapted from https://github.com/dib-lab/sourmash/blob/ca201cfc1824900bfb94aaadd94e163c70b7cf4e/src/sourmash/sig/__main__.py#L393

    Get hashes present in both sig1 and sig2, optionally using abundances from one of the signatures
    """
    hashes = set(sig1.minhash.hashes).intersection(sig2.minhash.hashes)

    if abundances_from is not None:
        hash_abundances = {hashval: sig1.minhash.hashes[hashval] for hashval in hashes}
    else:
        hash_abundances = dict.fromkeys(hashes)
    return hash_abundances


def get_kmers_in_seqfiles(
    seqfiles,
    query_hashvals,
    ksize,
    moltype,
    input_is_protein=True,
    seqout_fp=None,
    watermark=NOTIFY_EVERY_BP,
):

    n_seq = 0
    n = 0  # bp loaded
    m = 0  # bp in found sequences
    p = 0  # number of k-mers found

    found_kmers = []
    for filename in seqfiles:
        m, n = get_matching_hashes_in_file(
            filename,
            ksize,
            moltype,
            input_is_protein,
            query_hashvals,
            found_kmers,
            m,
            n,
            n_seq,
            seqout_fp,
            watermark,
        )
    return pd.DataFrame(
        found_kmers, columns=["kmer_in_seq", "kmer_in_alphabet", "hashval", "read_name"]
    )


def read_kmer_csv(csv, gene_name_tag="GN"):
    kmers = pd.read_csv(csv)
    kmers["gene_name"] = kmers["read_name"].str.extract(f"{gene_name_tag}:Z:([\w-]+)")
    return kmers


def celltype_cleaner(celltype):
    return celltype.lower().replace(" ", "_").replace("/", "-slash-")


def get_peptide_fasta(
    translate_base,
    channel,
    cell_barcode,
    is_aligned=True,
    double_aligned=True,
    channel_suffix=None,
):
    channel_suffix = "" if None else channel_suffix
    alignment_status = "aligned" if is_aligned else "unaligned"
    if double_aligned:
        alignment_status = "__".join([alignment_status, alignment_status])

    basename = (
        f"{channel}{channel_suffix}__{alignment_status}__{cell_barcode}__"
        + "coding_reads_peptides.fasta"
    )
    fasta = os.path.join(translate_base, basename)
    return fasta


def get_celltype_sig_path(celltype_sig_base, celltype_name, ksize, alphabet, scaled):
    """Per-celltype (merged individual signatures per cell type)

    Ends in .sig.sig for some reason
    """
    sketch_id = sig_utils.make_sketch_id(
        alpha=alphabet, ksize=ksize, style="scaled", value=scaled
    )
    sanitized_name = sig_utils.sanitize(celltype_name)
    sig_path = os.path.join(celltype_sig_base, sketch_id, f"{sanitized_name}.sig.sig")
    return sig_path


def get_cell_sig_path(
    cell_sig_base, cell_id, ksize, alphabet, scaled, add_ksize_to_sig_path=False
):
    """Per-cell-id signature. In a base directory for the sketch id"""
    sketch_id = sig_utils.make_sketch_id(
        alpha=alphabet, ksize=ksize, style="scaled", value=scaled
    )
    if add_ksize_to_sig_path:
        sig_path = os.path.join(cell_sig_base, sketch_id, str(ksize), f"{cell_id}.sig")
    else:
        sig_path = os.path.join(cell_sig_base, sketch_id, f"{cell_id}.sig")
    return sig_path


def get_diagnostic_kmers_for_cell(
    cell_id,
    cell_sig_base,
    cell_fasta_dir,
    predicted_celltype_name,
    #     ground_truth_celltype_name,
    celltype_sig_base,
    alphabet,
    ksize,
    scaled,
    double_aligned=False,
    verbose=False,
    gene_name_tag="GN",
    seqout_template=None,
    channel_suffix=None,
    add_ksize_to_sig_path=False,
):
    """

    seqout_template: str
        Where to write matching fasta files for kmers. Includes "{alignment_status}"
        in the name to put aligned/unaligned
    """
    celltype_sig_path = get_celltype_sig_path(
        celltype_sig_base,
        predicted_celltype_name,
        ksize=ksize,
        alphabet=alphabet,
        scaled=scaled,
    )
    cell_sig_path = get_cell_sig_path(
        cell_sig_base,
        cell_id,
        ksize=ksize,
        alphabet=alphabet,
        scaled=scaled,
        add_ksize_to_sig_path=add_ksize_to_sig_path,
    )

    channel, cell_barcode = cell_id.split("__")
    aligned_fasta = get_peptide_fasta(
        cell_fasta_dir,
        channel,
        cell_barcode,
        is_aligned=True,
        double_aligned=double_aligned,
        channel_suffix=channel_suffix,
    )
    unaligned_fasta = get_peptide_fasta(
        cell_fasta_dir,
        channel,
        cell_barcode,
        is_aligned=False,
        double_aligned=double_aligned,
        channel_suffix=channel_suffix,
    )

    fastas = {"aligned": aligned_fasta, "unaligned": unaligned_fasta}

    cell_sigobj = sourmash.load_one_signature(
        cell_sig_path, ksize=ksize, select_moltype=alphabet
    )
    celltype_sigobj = sourmash.load_one_signature(
        celltype_sig_path, ksize=ksize, select_moltype=alphabet
    )

    if verbose:
        overlap_stats = overlap(
            cell_sigobj, celltype_sigobj, moltype=alphabet, ksize=ksize
        )
        print(overlap_stats)

    intersecting_hashes = get_intersecting_hashes(
        cell_sigobj, celltype_sigobj, abundances_from=cell_sigobj
    )
    dfs = []
    for alignment_status, fasta in fastas.items():
        if seqout_template is not None:
            seqout_filename = seqout_template.format(alignment_status=alignment_status)
            dirname = os.path.dirname(seqout_filename)
            if not os.path.exists(dirname):
                os.makedirs(dirname)

            seqout_fp = open(seqout_filename, "wt")
        else:
            seqout_fp = None

        try:
            df = get_kmers_in_seqfiles(
                [fasta],
                query_hashvals=intersecting_hashes.keys(),
                ksize=ksize,
                moltype=alphabet,
                seqout_fp=seqout_fp,
            )
        except FileNotFoundError:
            logging.error(f"Could not open {fasta}")
            return

        if seqout_template is not None:
            seqout_fp.close()

        df["alignment_status"] = alignment_status
        dfs.append(df)
    kmers = pd.concat(dfs)
    kmers["gene_name"] = kmers["read_name"].str.extract(f"{gene_name_tag}:Z:([\w-]+)")
    kmers["abundances"] = kmers.hashval.map(intersecting_hashes)

    kmers["cell_id"] = cell_id
    kmers["predicted_celltype_name"] = predicted_celltype_name
    return kmers


if __name__ == "__main__":
    sys.exit(main())
