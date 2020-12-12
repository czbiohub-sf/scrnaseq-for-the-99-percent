
import re
import os
import glob
import copy
import logging


from joblib import Parallel, delayed
import numpy as np
import sourmash
from sourmash import sourmash_args
import pandas as pd
from sourmash.sourmash_args import FileOutput

import sig_utils

logging.basicConfig(filename='create_sourmash_commands.log', level=logging.INFO)


def get_cell_id_from_fasta(fasta_basename_noextension):
    # Remove coding reads
    name = re.sub(r"__noncoding_reads_nucleotides|__coding_reads_nucleotides", "", fasta_basename_noextension)
    sub1 = re.sub(
        r"__aligned__aligned__|__unaligned__unaligned__", "__", name
    )
    cell_id = re.sub(
        r"__aligned__|__unaligned__", "__", sub1
    ).replace("__coding_reads_peptides", "")
    return cell_id
    


def make_sourmash_compute_commands(
    fastas: str, 
    compute_txt:str, 
    sig_output_folder:str, 
    ksizes=sig_utils.KSIZES, 
    scaled=10, 
    n_jobs=96,
    input_is_protein=True,
    cell_id_as_name=True
):
    """
    PARAMS:
    -------
    fastas: fastas you want to make sigs from 
    compute_txt: name of file to write commands to
    sig_output_folder: output folder for new sigs
    cell_id_as_name: use the cell ID as the name. Not useful for predictorthologs because need aligned/unaligned in the signature name
    """
    
    ksizes_flag = f'--ksizes {",".join(map(str, ksizes))}'
    with open(compute_txt, "w") as fh:
        for fasta in fastas:
            f = fasta
            f_base = os.path.basename(f)
            f_base_noext = os.path.splitext(f_base)[0]
            cell_id = get_cell_id_from_fasta(f_base_noext)
            
            if cell_id_as_name:
                name = cell_id
            else:
                name = f_base_noext

            try:
                output = os.path.join(
                    sig_output_folder,
                    f_base_noext + ".sig"
                )
                # removed --protein from command since we only use dayhoff signatures, and it will make loading the signatures ~2x faster (I think)
                if input_is_protein:
                    moltype_flag = '--input-is-protein --protein --dayhoff --hp'
                else: 
                    moltype_flag = '--dna'
                command = f"sourmash compute --quiet --track-abundance --scaled {scaled} {ksizes_flag} {moltype_flag} {fasta} --name '{name}' -o {output} \n"
                fh.write(command)

            except KeyError:
                pass
    print(f"parallel --progress --eta --jobs {n_jobs} < {compute_txt}")

# use sig_utils to merge aligned and unaligned sigs
    

def remove_ribosomal(
    ksize: int, 
    cell_id: str,
    original_sig_path: str,
    ribosomal_sig_path: str, 
    remove_ribosomal_output_dir: str, 
    moltype="dayhoff",
    force=False
    
):
    """
    * probably should do this before building the SBTs
    """
    sketch_id = sig_utils.make_sketch_id(
        moltype, ksize, style='scaled', value=10
    ) 
    remove_ribosomal_sketch_id_dir = os.path.join(
        remove_ribosomal_output_dir, 
        sketch_id
    ) 
    output_sig = os.path.join(
        remove_ribosomal_sketch_id_dir, 
        os.path.basename(original_sig_path)
    )
    
    if not force and os.path.exists(output_sig):
        return    

    
    os.makedirs(remove_ribosomal_sketch_id_dir, exist_ok=True)
    
    try:
        original_sig = sourmash.load_one_signature(
            original_sig_path,
            ksize=ksize,
            select_moltype=moltype
        )
    except ValueError:
        logging.error(f"Could not load {original_sig_path} with ksize={ksize} and moltype={moltype}")
        return
        
    original_sig_copy = copy.deepcopy(original_sig)
    flattened_sig = sig_utils._flatten(
        original_sig_copy, 
        ksize, 
        moltype,
        cell_id,
        md5=None, 
        name=None
    )
    
    subtracted_sig = sig_utils._subtract(
        flattened_sig, 
        ribosomal_sig_path, 
        ksize,
        cell_id,
        moltype=moltype
    )
    
    intersected_sig = sig_utils._intersect(
        ksize,  
        [subtracted_sig, original_sig], 
        original_sig,
        cell_id, # to rename 
        moltype=moltype,
    )
       
    try:
        with FileOutput(output_sig, 'wt') as fp:
            sourmash.save_signatures([intersected_sig], fp=fp)
            logging.info(f"intersected sig successfully save {intersected_sig} for cell_id {cell_id}")
    except: 
        logging.debug(f"fileoutput failed for {original_sig_path}, {ksize}, {cell_id}")
        
    
def make_sourmash_sbts(
    sbt_base_dir: str, 
    ksizes: list, 
    merged_sigs_output_folder: str, 
    index_command_file_name: str, 
    moltype='dayhoff',
    sbt_prefix='mouse',
    n_jobs=22,
    
):
    """
    sbt_base_dir: this is where the output SBTs will be
    ksizes: k size: (will build an SBT for each ksize provided)
    merged_sigs_output_folder: Folder of signatures to build the SBT off of 
    index_command_file_name: name of textfile to output command to 
    """
    with open(index_command_file_name, "w") as fh:
        for k in ksizes:
            sketch_id = sig_utils.make_sketch_id(moltype, k, style='scaled', value=10) 
            if moltype.lower() == 'dna':
                moltype_args = f'--{moltype.lower()} --no-protein --no-dayhoff'
            else:
                moltype_args = f'--{moltype} --no-dna'
            flags = f"--traverse-directory --quiet {moltype_args}"
            command = f"sourmash index -k {k} {flags} {sbt_base_dir}/{sbt_prefix}__alphabet-{moltype}__ksize-{k}.sbt.zip {merged_sigs_output_folder}/{sketch_id} \n"
            fh.write(command)

    print(f"parallel --progress --eta --jobs {n_jobs} < {index_command_file_name}")
    return index_command_file_name


def make_sourmash_search_commands(
    search_dir, 
    merged_sigs_dir, 
    sbt_base_dir, 
    k_sizes, 
    scaled_sizes=[10,],
    cell_ids=[],
    sbt_template_basename=r"mouse_ksize_{k}.sbt.zip", 
    query_sig_files=False, 
    containment=True,
    moltype="dayhoff",
    num_results=3,
    n_jobs=96
):    
    """
    PARAMS:
    -------
    search_dir : where you'll write search results csv file
    merged_sigs_dir: merged sigs that you are searching agaisnt SBTs
    sbt_base_dir: this will be the directory to your sbts
    k_sizes: used for grouping signatures/SBTs to search
    scaled_sizes=[10,] : inputs 
    cell_ids: if interested in only certain cell_ids in groupings
    sbt_template_basename=r"mouse_ksize_{k}.sbt.zip": what template we use to name the SBTs 
    query_sig_files: this is the signature files to be queryed
    containment=True: if you want to search with containment instead of similarity
    moltype=dayhoff: molecule or alphabet type. either 'dayhoff', 'protein', or 'dna'
    num_results=5: Number of search results to return. Fewer means search goes faster
    """
#     moltype_lower = moltype.lower()
    valid_moltypes = 'dayhoff', 'protein', 'dna'
    if not isinstance(moltype, str) or moltype.lower() not in valid_moltypes:
        raise ValueError(f"{moltype} is not a valid molecule type. Only one of {valid_moltypes} are allowed")
    
    sourmash_search_command_file = os.path.join(search_dir, f"sourmash_search_commands_{moltype}.txt")
    with open(sourmash_search_command_file, "w") as fh:
        for k in k_sizes:
            for scale in scaled_sizes:
                sketch_id = sig_utils.make_sketch_id(moltype, k, style='scaled', value=scale)
                sketch_id_merged_sig_dir = os.path.join(merged_sigs_dir, sketch_id)
                sketch_id_search_dir = os.path.join(search_dir, sketch_id) # where we'll write search results to

                if not os.path.exists(sketch_id_search_dir):
                    os.makedirs(sketch_id_search_dir)

                if not cell_ids:
                    query_sig_files = [
                        os.path.join(sketch_id_merged_sig_dir, sig) for sig in 
                        glob.glob(os.path.join(sketch_id_merged_sig_dir, "*.sig"))
                    ]
                else:
                    # if cell id provided then create query sig files 
                    query_sig_files = [
                        os.path.join(sketch_id_merged_sig_dir, f"{cell_id}.sig") 
                        for cell_id in cell_ids
                    ]
                    
                sbt_index = os.path.join(sbt_base_dir, sbt_template_basename.format(k=k, moltype=moltype))
                for sig in query_sig_files:
                    #complete_sig_path = os.path.join(sketch_id_merged_sig_dir, sig)
                    basename = os.path.basename(sig)
                    output_csv = os.path.join(sketch_id_search_dir, basename.replace(".sig", ".csv"))
                    if moltype.lower() == 'dna':
                        moltype_args = f'--{moltype} --no-protein --no-dayhoff'
                    else:
                        moltype_args = f'--{moltype} --no-dna'
                    flags = f"--quiet {moltype_args} --num-results {num_results} --threshold 0.0001 -k {k}"
                    if containment:
                        flags += " --containment"
                    command = f"sourmash search {flags} --output {output_csv} {sig} {sbt_index}\n"
                    fh.write(command)
                                
    print(f"parallel --progress --eta --jobs {n_jobs} < {sourmash_search_command_file}")
    return sourmash_search_command_file

    
def get_cells_with_both_aligned_unaligned_sigs_df(sig_folder):
    globber = os.path.join(sig_folder, "*.sig")
    
    df = pd.Series(glob.glob(globber), name='fullpath').to_frame()
    df['basename'] = df['fullpath'].map(os.path.basename)
    df["cell_id"] = df.basename.str.replace(
        r"__aligned__aligned__|__unaligned__unaligned__", "__"
    ).str.replace(
        '__aligned__|__unaligned', '__').str.split(
        "__coding_reads_peptides").str[0]
    #df.head()

    df['fasta_id'] = df.basename.str.split('__coding_reads').str[0]
    df['channel'] = df.cell_id.str.split('__').str[0]
    df['cell_barcode'] = df.cell_id.str.split('__').str[1]
    print(df.shape)

    both_aligned_unaligned = df.groupby("cell_id").filter(lambda x: len(x)==2)
    print(both_aligned_unaligned.shape)
    both_aligned_unaligned.sort_values("cell_id").head()
    return both_aligned_unaligned



def make_merged_sigs_df(
    merged_sigs_dir,
    globber_template='{merged_sigs_dir}/*/*.sig'
):
    globber = globber_template.format(merged_sigs_dir=merged_sigs_dir)
    merged_sigs = df = pd.Series(glob.glob(globber), name='fullpath'
    ).to_frame()
    merged_sigs['basename'] = merged_sigs['fullpath'].map(os.path.basename)
    merged_sigs["cell_id"] = merged_sigs.basename.str.split(".sig").str[0]
    merged_sigs['channel'] = merged_sigs.cell_id.str.split('__').str[0]
    merged_sigs['cell_barcode'] = merged_sigs.cell_id.str.split('__').str[-1]
    merged_sigs["alphabet"] = merged_sigs.fullpath.str.split("alphabet-").str[-1].str.split("__").str[0]
    merged_sigs["ksize"] = merged_sigs.fullpath.str.split("ksize-").str[-1].str.split("__").str[0].astype(int)
    merged_sigs["scaled"] = merged_sigs.fullpath.str.split("scaled-").str[-1].str.split("/").str[0].astype(int)
    merged_sigs['sketch_id'] = merged_sigs.fullpath.str.split('/').str[-2]
    return merged_sigs


def join_sigs_with_ontologies(
    merged_sigs,
    metadata, # where to join cell_onotology class from, ex: one2one.obs
    metadata_cell_ontology_cols= ['cell_ontology_class', "broad_group"],
    metadata_join_cols = ["cell_barcode", "channel"],
    sample_from_col = "broad_group",
    cell_ontology_groups = [
        'Dendritic', 
        'Macrophage',
        'Monocyte',
        'B cell',
        'T cell',
        'Natural Killer T cell',
        'Natural Killer',
        'Alveolar Epithelial Type 2',
        'Ciliated',
        'Fibroblast',
        'Smooth Muscle and Myofibroblast'
    ],
):
    TOTAL_N_SKETCH_IDS = merged_sigs.sketch_id.nunique()

#     import pdb; pdb.set_trace()
    # To do apples-to-apples comparisons, filter for only cell IDs present in all sketches
    merged_sigs = merged_sigs.groupby('cell_id').filter(lambda x: x.sketch_id.nunique() == TOTAL_N_SKETCH_IDS)
    
    # merge on ontology:
    merged_sigs_w_ontology = merged_sigs.merge(
        metadata[metadata_cell_ontology_cols + metadata_join_cols], 
        left_on=metadata_join_cols,
        right_on=metadata_join_cols
    )
    
    # subset to broad tissue types
    if cell_ontology_groups:
        merged_sigs_w_ontology_subset = merged_sigs_w_ontology[
            merged_sigs_w_ontology[sample_from_col].isin(cell_ontology_groups)
        ]
        return merged_sigs_w_ontology_subset
    
    else:
        print("not taking subsets")
        return merged_sigs_w_ontology.dropna(subset=metadata_cell_ontology_cols)
    
    
def subsample_sig_df_ontologies(
    merged_sigs_w_ontology_subset,
    n_samples=10,
    sample_from_col="broad_group",
):
        
    # This is a pretty "dumb way" to do this but it's what I can think of right now.. will do something more elegant later
#     Groupby sketch id, but only do this once since each cell should be represented once per sketch id
    one_sketch_id = merged_sigs_w_ontology_subset.head(1)['sketch_id'].values[0]
    df = merged_sigs_w_ontology_subset.query('sketch_id == @one_sketch_id')
    sampled_number_of_cell_ids = df.groupby(sample_from_col, observed=True).apply(
        lambda x: x.sample(n_samples, random_state=0, replace=False) if len(x) > n_samples else x
    )
    
    merged_sigs_w_ontology_subset = merged_sigs_w_ontology_subset.query("cell_id in @sampled_number_of_cell_ids.cell_id.values")
    
    return merged_sigs_w_ontology_subset


def sample_sigs_from_ontologies(
    merged_sigs_dir,
    metadata, # where to join cell_onotology class from, ex: one2one.obs
    metadata_cell_ontology_cols= ['cell_ontology_class', "broad_group", 'compartment_group'],
    metadata_join_cols = ["cell_barcode", "channel"],
    sample_from_col = "broad_group",
    cell_ontology_groups = [
        'Dendritic', 
        'Macrophage',
        'Monocyte',
        'B cell',
        'T cell',
        'Natural Killer T cell',
        'Natural Killer',
        'Alveolar Epithelial Type 2',
        'Ciliated',
        'Fibroblast',
        'Smooth Muscle and Myofibroblast'
    ],
    n_samples=10,
):

    merged_sigs = make_merged_sigs_df(merged_sigs_dir)
    
    mouse_sigs_ontologies = join_sigs_with_ontologies(
        merged_sigs, 
        metadata=metadata,
        metadata_cell_ontology_cols=metadata_cell_ontology_cols,
        metadata_join_cols=metadata_join_cols,
        sample_from_col=sample_from_col,
        cell_ontology_groups=cell_ontology_groups
    )
    
    mouse_sigs_ontologies_sampled = subsample_sig_df_ontologies(
        mouse_sigs_ontologies, 
        n_samples=n_samples,
        sample_from_col=sample_from_col,
    )
    return mouse_sigs_ontologies_sampled