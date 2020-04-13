DIR=/home/olga/data_sm/tabula-microcebus/from_stanford_sherlock/oak/stanford/groups/krasnow/ktrav/HLCA/data10x/sequencing_runs/sequencing_runs/171205_A00111_0088_BH55NYDMXX/alignments
WORK=${HOME}/data_sm/nextflow-intermediates-ndnd/
OUT=/mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--hlca/
PEPTIDE_FASTA=/home/olga/data_lg/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_2759.fasta.gz


SOFTLINK_DIR=/home/olga/data_sm/tabula-microcebus/rawdata/tenx/lung-soft-links//human

COMMON_FLAGS=		nf-core/kmermaid \
		-r olgabot/10x-tgz \
		-latest \
		-profile docker \
		-work-dir ${WORK} \
		-resume \
		--tenx_min_umi_per_cell 5001 \
		-with-tower \
		--track_abundance \
		--peptide_fasta ${PEPTIDE_FASTA} \
		--extract_coding_peptide_ksize 9 \
		--extract_coding_jaccard_threshold 0.95 \
		--molecules dna,protein,dayhoff \
		--ksizes 21,24,27,33,36,39,42,45,51 \
		--log2_sketch_sizes 14 \


firstrun:
	nextflow run \
		${COMMON_FLAGS} \
		--tenx_tgz ${DIR}/'**.tgz' \
		--outdir ${OUT}/$@


normal_lung_only:
	nextflow run \
		${COMMON_FLAGS} \
		--tenx_tgz ${SOFTLINK_DIR}/'**.tgz' \
		--outdir ${OUT}/$@
