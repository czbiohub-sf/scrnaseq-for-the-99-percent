
#! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

Cribbed from https://github.com/dib-lab/sourmash/pull/724/
"""
import sys
import argparse
import sourmash
from sourmash import MinHash
from sourmash import sourmash_args
from sourmash._minhash import hash_murmur
import screed
import csv
from sourmash.logging import notify, error
from sourmash.cli.utils import add_construct_moltype_args, add_ksize_arg
from sourmash.sourmash_args import calculate_moltype
from sencha.sequence_encodings import encode_peptide, AMINO_ACID_SINGLE_LETTERS



NOTIFY_EVERY_BP = 1e7


def get_kmer_moltype(sequence, start, ksize, moltype, input_is_protein):
    kmer_in_seq = sequence[start : start + ksize]
    if moltype == "DNA":
        # Get reverse complement
        kmer_rc = screed.rc(kmer)
        if kmer_in_seq > kmer_rc:  # choose fwd or rc
            kmer_encoded = kmer_rc
        else:
            kmer_encoded = kmer_in_seq
    elif input_is_protein:
        kmer_encoded = encode_peptide(kmer_in_seq, moltype)
    elif not input_is_protein:
        raise NotImplementedError(
            "Currently cannot translate DNA to protein sequence"
        )
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

        kmer_encoded, kmer_in_seq = get_kmer_moltype(sequence, start, ksize, moltype, input_is_protein)

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
            found_kmers.append([kmer_in_seq, kmer_encoded, hashval, record['name']])

            # write out sequence
            if seqout_fp:
                seqout_fp.write(">{}|hashval:{}|kmer:{}|kmer_encoded:{}\n{}\n".format(
                    record.name, hashval, kmer_in_seq, kmer_encoded, record.sequence))
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
    add_ksize_arg(p)
    add_construct_moltype_args(p)
    args = p.parse_args()

    # set up the outputs.
    seqout_fp = None
    if args.output_sequences:
        seqout_fp = open(args.output_sequences, "wt")

    kmerout_fp = None
    if args.output_kmers:
        kmerout_fp = open(args.output_kmers, "wt")
        kmerout_w = csv.writer(kmerout_fp)
        kmerout_w.writerow(["kmer_in_sequence", "kmer_in_alphabet", "hashval", "read_name"])

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
    sigobj = sourmash.load_one_signature(args.query, ksize=args.ksize, select_moltype=moltype)
    query_hashvals = set(sigobj.minhash.get_mins())
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


if __name__ == "__main__":
    sys.exit(main())
