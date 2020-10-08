TGZ_DIR=/data_sm/olga/tabula-microcebus/rawdata/tenx/lung-soft-links/human
TRANSLATE_FASTA=/data_lg/olga/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz
OUTDIR_BASE=/data_sm/olga/tabula-microcebus/analyses/kmermaid/tenx-tgz--hlca-v2

normal-lung:
	nextflow run \
                nf-core/kmermaid \
		-r pv-patch-2 \
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
		-work-dir ~/data_sm/nextflow-work-lrrr/ \
		 --tenx_min_umi_per_cell 1000 \
		-with-tower \
		-resume
