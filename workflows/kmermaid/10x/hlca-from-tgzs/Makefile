TGZ_DIR=/mnt/data_sm/olga/tabula-microcebus/rawdata/tenx/lung-soft-links/human
TRANSLATE_FASTA=/mnt/data_lg/olga/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz
OUTDIR_BASE=/mnt/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--hlca-v2

normal-lung:
	/home/pranathi/nextflow run \
                nf-core/kmermaid \
		-r pranathivemuri-patch-1 \
		-latest \
		--tenx_tgz ${TGZ_DIR}'/P1**.tgz' \
		--skip_compute \
		--skip_compare \
		--molecules dayhoff \
		--ksizes 27,30,33,51,54,57,60 \
		--sketch_scaled 10 \
		--track_abundance \
		--reference_proteome_fasta ${TRANSLATE_FASTA} \
		--outdir ${OUTDIR_BASE}--$@ \
		-work-dir /mnt/data_lg/pranathi/nextflow-work-lrrr/ \
		 --tenx_min_umi_per_cell 1000 \
		-with-tower \
		-resume --skip_multiqc

normal-lung-remove-riborna:        
	/home/pranathi/nextflow run \
                nf-core/kmermaid \
                -r pranathivemuri-patch-1 \
                -latest \
                --tenx_tgz ${TGZ_DIR}'/P1**.tgz' \
                --skip_compute \
                --skip_compare \
                --molecules dayhoff \
                --ksizes 27,30,33,51,54,57,60 \
                --sketch_scaled 10 \
                --track_abundance \
                --reference_proteome_fasta ${TRANSLATE_FASTA} \
                --outdir ${OUTDIR_BASE}--$@ \
                -work-dir /mnt/data_lg/pranathi/nextflow-work-lrrr/ \
                 --tenx_min_umi_per_cell 1000 \
                -with-tower --remove_ribo_rna --save_non_rrna_reads \
                -resume --skip_multiqc


normal-lung-remove-riborna-p3:
	/home/pranathi/nextflow run \
                nf-core/kmermaid \
                -r pranathivemuri-patch-1 \
                -latest \
                --tenx_tgz ${TGZ_DIR}'/P3**.tgz' \
                --skip_compute \
                --skip_compare \
                --molecules dayhoff \
                --ksizes 27,30,33,51,54,57,60 \
                --sketch_scaled 10 \
                --track_abundance \
                --reference_proteome_fasta ${TRANSLATE_FASTA} \
                --outdir ${OUTDIR_BASE}--$@ \
                -work-dir /mnt/data_lg/pranathi/nextflow-work-lrrr/ \
                 --tenx_min_umi_per_cell 1000 \
                -with-tower --remove_ribo_rna --save_non_rrna_reads \
                -resume --skip_multiqc
