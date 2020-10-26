TRANS_PROTEOME=/home/olga/data_lg/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz
TGZS_DIR=/mnt/data_sm/olga/immune-evolution/rawdata/mouse/tgzs/*LUNG*
OUTDIR_BASE=/mnt/data_sm/olga/immune-evolution/pipeline-results/mouse/kmermaid/lung
WORK_DIR=/mnt/data_lg/lekha/mouse

mouse:
	nextflow run nf-core/kmermaid \
                -r dev \
                -latest \
                -profile docker \
                -resume \
                -with-tower \
                --molecules dayhoff \
		--skip_compare \
		--skip_compute \
		--track_abundance \
                --ksizes 27,30,33,36,39,42,45,48,51,54,57,60,63 \
                --sketch_scaled 10 \
                -work-dir ${WORK_DIR} \
                --translate_peptide_ksize 11 \
                --reference_proteome_fasta ${TRANS_PROTEOME} \
                        --tenx_tgz ${TGZS_DIR}'/**.tgz' \
			--tenx_min_umi_per_cell 1000 \
		--outdir ${OUTDIR_BASE}--$@


mouse--remove-ribo:
	
	nextflow run nf-core/kmermaid \
		-r dev \
		-latest \
		-profile docker \
		-resume \
		-with-tower \
		--molecules dayhoff \
		--skip_compare \
		--skip_compute \
		--track_abundance \
		--ksizes 27,30,33,36,39,42,45,48,51,54,57,60,63 \
		--sketch_scaled 10 \
		--remove_ribo_rna \
		--save_non_rrna_reads \
		-work-dir ${WORK_DIR} \
		--translate_peptide_ksize 8 \
		--translate_jaccard_threshold 0.05 \
		--reference_proteome_fasta ${TRANS_PROTEOME} \
			--tenx_tgz ${TGZS_DIR}'/**.tgz' \
	        	--tenx_min_umi_per_cell 1000 \
		--outdir ${OUTDIR_BASE}--$@
