REFERENCE_PROTEOME_FASTA=/home/olga/data_lg/czbiohub-reference/uniprot/taxonomy_Mammalia_9MAMM_40674_plus_HLrhiFer5_Bat1K.fasta.gz
#BAM_DIR=/mnt/ibm_sm/home/phoenix/batlas/BatBams/immune-tissue-softlinks
READ_SINGLES='/home/olga/data_sm/immune-evolution/pipeline-results/bat/kmermaid/run-1/10x-fastqs/per-cell/*LU*/*.fastq.gz'
OUTDIR_BASE=/mnt/ibm_sm/olga/immune-evolution/pipeline-results/bat/kmermaid/
WORK_DIR=/mnt/ibm_sm/home/phoenix/nextflow-work

COMMON_FLAGS= \
		-latest \
		-resume \
		-with-tower \
		--molecules protein,dayhoff \
		--sketch_scaled 10 \
		--reference_proteome_fasta ${REFERENCE_PROTEOME_FASTA} \
		--tenx_min_umi_per_cell 1000 \
		--read_singles ${READ_SINGLES} \
		--track_abundance \
		--skip_compute \
		--skip_compare \
		--skip_multiqc \
#		--ksizes 27,30,33,36,39,42,45,48,51,54,57,60,63

run-1:
	nextflow run nf-core/kmermaid \
		-r phoenix/no-multiqc \
		-work-dir ${WORK_DIR} \
		${COMMON_FLAGS} \
		--translate_peptide_ksize 11 \
		--outdir ${OUTDIR_BASE}/$@


ksize_8:
	nextflow run nf-core/kmermaid \
		-r phoenix/no-multiqc \
		-work-dir ${WORK_DIR} \
		${COMMON_FLAGS} \
		--translate_peptide_ksize 8 \
		--translate_jaccard_threshold 0.05 \
		--remove_ribo_rna \
		--save_non_rrna_reads \
		--outdir ${OUTDIR_BASE}/$@
