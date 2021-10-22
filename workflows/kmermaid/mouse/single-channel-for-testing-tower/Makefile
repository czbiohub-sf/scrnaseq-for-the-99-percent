

all:
	nextflow run nf-core/kmermaid \
		-r olgabot-patch-3 \
		-profile docker \
		--reference_proteome_fasta s3://immune-evolution/databases/uniprot/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz \
		--bam 's3://czb-tabula-muris-senis/10x/18_month/MACA_18m_F_LUNG_50/*.bam' \
		--skip_compute \
		--skip_compare \
		--sketch_scaled 10 \
		--tenx_min_umi_per_cell 1000 
