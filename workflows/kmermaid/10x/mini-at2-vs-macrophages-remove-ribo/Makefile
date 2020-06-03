DATA_DIR=/mnt/data_sm/olga/tabula-microcebus/rawdata/tenx/mini-test-at2-vs-macrophages-soft-links
WORK=${HOME}/data_sm/nextflow-intermediates-lrrr/
OUT=/mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/mini-test-at2-vs-macrophages/
PEPTIDE_FASTA=/home/olga/data_lg/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_2759.fasta.gz
TRANSLATE_FASTA="/home/olga/data_lg/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz"


dev:
	nextflow run \
		nf-core/kmermaid \
		-r dev \
		-latest \
		-profile docker \
		--tenx_tgz ${DATA_DIR}/'**.tgz' \
		-work-dir ${WORK} \
		--outdir ${OUT}/$@ \
		-resume \
		--tenx_min_umi_per_cell 5001 \
		-with-tower \
		--track_abundance \
		--peptide_fasta ${PEPTIDE_FASTA} \
		--translate_peptide_ksize 9 \
		--translate_jaccard_threshold 0.95 \
		--molecules dna,protein,dayhoff,hp \
		--ksizes 21,24,27,33,36,39,42,45,51,57,63,69,75,81,87,93,99,105,111,117 \
		--log2_sketch_sizes 10



remove_ribo:
	nextflow run \
		nf-core/kmermaid \
		-r olgabot/remove-ribosomal-reads \
		-latest \
		-profile docker \
		--tenx_tgz ${DATA_DIR}/'**.tgz' \
		-work-dir ${WORK} \
		--outdir ${OUT}/$@ \
		-resume \
		--tenx_min_umi_per_cell 5001 \
		-with-tower \
		--track_abundance \
		--peptide_fasta ${PEPTIDE_FASTA} \
		--translate_peptide_ksize 9 \
		--translate_jaccard_threshold 0.95 \
		--molecules dna,protein,dayhoff,hp \
		--ksizes 21,24,27,33,36,39,42,45,51,57,63,69,75,81,87,93,99,105,111,117 \
		--log2_sketch_sizes 10 \
		--remove_ribo_rna \
		--save_non_rrna_reads


remove_ribo_ksize_mammals:
	nextflow run \
		nf-core/kmermaid \
		-r olgabot/propagate-peptide-ksize \
		-latest \
		-profile docker \
		--tenx_tgz ${DATA_DIR}/'**.tgz' \
		-work-dir ${WORK} \
		--outdir ${OUT}/$@ \
		-resume \
		--tenx_min_umi_per_cell 5001 \
		-with-tower \
		--track_abundance \
		--peptide_fasta ${TRANSLATE_FASTA} \
		--translate_peptide_ksize 9 \
		--translate_jaccard_threshold 0.95 \
		--molecules dna,protein,dayhoff \
		--ksizes 51 \
		--log2_sketch_sizes 10 \
		--remove_ribo_rna \
		--save_non_rrna_reads

