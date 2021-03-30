import glob
import os
import shutil

from tqdm import tqdm

NUCLEOTIDE_SKETCH_IDS = ["alphabet-DNA__ksize-21__scaled-10"]

PEPTIDE_SKETCH_IDS = [
    "alphabet-protein__ksize-30__scaled-10",
    "alphabet-dayhoff__ksize-51__scaled-10",
]

SHARED_CELLTYPES = [
    "Alveolar Epithelial Type 2",
    "B cell",
    "Capillary",
    "Dendritic",
    "Fibroblast",
    "Macrophage",
    "Monocyte",
    "Natural Killer T cell",
    "Smooth Muscle and Myofibroblast",
    "T cell",
]


def copy_fastas(
    fasta_outdir_base,
    input_fasta_dir,
    select_cell_ids=None,
    dryrun=False,
    molecule_types=("peptides", "nucleotides"),
    upper=False,
    cell_id_fun=None
):
    """Copy k-mer signatures under the same folder per sketch id
    
    
    cell_id_fun : function
        String manipulation function to clean cell id from the basename of the 
        sketch to have it match the list of select cell ids
    """
    for moltype in molecule_types:
        print(f"Copying {moltype}")
        outdir = os.path.join(fasta_outdir_base, moltype)
        if not dryrun and not os.path.exists(outdir):
            os.makedirs(outdir)

        fasta_glob = os.path.join(input_fasta_dir, f"*_coding_reads_{moltype}.fasta")

        for fasta_original in tqdm(glob.iglob(fasta_glob)):

            basename = os.path.basename(fasta_original)
            cell_id = (
                basename.replace("__aligned__aligned", "")
                .replace("__aligned__", "__")
                .replace("__unaligned__unaligned", "")
                .replace("__unaligned__", "__")
                .split("__coding")[0]
            )
            if upper:
                cell_id = cell_id.upper()
            
            # Apply cell id cleaning function if specified    
            if cell_id_fun is not None:
                cell_id = cell_id_fun(cell_id)
#             import pdb; pdb.set_trace()
                
            if select_cell_ids is not None and cell_id not in select_cell_ids:
                continue

            fasta_newplace = os.path.join(outdir, basename)
            if dryrun:
                print(f"Copy:\n{fasta_original}\n--> {fasta_newplace}")
            if not dryrun and not os.path.exists(fasta_newplace):
                shutil.copy(fasta_original, fasta_newplace)


def _per_sketch_id_copy_sketches(
    sketch_ids,
    input_sketch_dir,
    pre_sketch_id_outdir,
    select_cell_ids=None,
    dryrun=False,
    cell_id_fun=None
):
    """Copy k-mer signatures under the same folder per sketch id
    
    cell_id_fun : function
        String manipulation function to clean cell id from the basename of the 
        sketch to have it match the list of select cell ids
    """
    for sketch_id in sketch_ids:
        print(f"Copying {sketch_id}")
        outdir = os.path.join(pre_sketch_id_outdir, sketch_id)
        if not dryrun and not os.path.exists(outdir):
            os.makedirs(outdir)

        sigfile_glob = os.path.join(input_sketch_dir, sketch_id, "*.sig")

        for sigfile_original in tqdm(glob.iglob(sigfile_glob)):
            basename = os.path.basename(sigfile_original)
            
            # Don't need to check for aligned/unaligned here because the signatures are already merged
            cell_id = basename.split(".")[0]
            if cell_id_fun is not None:
                cell_id = cell_id_fun(cell_id)
#             import pdb; pdb.set_trace()
            
            if select_cell_ids is not None and cell_id not in select_cell_ids:
                continue

            # Set to the new, cleaned cell id if applicable
            sigfile_newplace = os.path.join(outdir, f"{cell_id}.sig")
            if dryrun:
                print(f"Copy:\n{sigfile_original}\n--> {sigfile_newplace}")
            if not dryrun and not os.path.exists(sigfile_newplace):
                shutil.copy(sigfile_original, sigfile_newplace)


def copy_nucleotide_peptide_sketches(
    peptide_sketch_dir,
    nucleotide_sketch_dir,
    pre_sketch_id_outdir,
    nucleotide_sketch_ids=NUCLEOTIDE_SKETCH_IDS,
    peptide_sketch_ids=PEPTIDE_SKETCH_IDS,
    select_cell_ids=None,
    dryrun=False,
    cell_id_fun=None
):
    """Copy both nucleotide and peptide sketches per sketch id
    
    
    cell_id_fun : function
        String manipulation function to clean cell id from the basename of the 
        sketch to have it match the list of select cell ids
dd    """
    if nucleotide_sketch_dir is not None:
        _per_sketch_id_copy_sketches(
            nucleotide_sketch_ids,
            nucleotide_sketch_dir,
            pre_sketch_id_outdir,
            select_cell_ids,
            dryrun,
            cell_id_fun=cell_id_fun
        )
    if peptide_sketch_dir is not None:
        _per_sketch_id_copy_sketches(
            peptide_sketch_ids,
            peptide_sketch_dir,
            pre_sketch_id_outdir,
            select_cell_ids,
            dryrun,
            cell_id_fun=cell_id_fun
        )
