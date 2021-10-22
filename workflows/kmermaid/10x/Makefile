


run:
	nextflow run \
                nf-core/kmermaid \
                -r olgabot/khtools-extract-coding \
		-latest \
		-profile czbiohub_local,docker \
		--bam  ~/data_sm/tabula-microcebus/rawdata/aligned_micmur3/10x/antoine__180917_A00111_0211_AHGKCVDMXX/ANTOINE_LUNG_P3/ANTOINE_LUNG_P3/outs/possorted_genome_bam.bam \
		--custom_config_base ~/code/nf-core/configs \
		--molecules dna,protein,dayhoff \
		--peptide_fasta ~/data_lg/czbiohub-reference/ensembl/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz \
		--ksizes 21,24,27,33,36 \
		--log2_sketch_sizes 10,12,14,16 \
		--bloomfilter_tablesize '1e8' \
		--extract_coding_peptide_ksize 11 \
		--extract_coding_peptide_molecule dayhoff \
		--outdir ~/data_sm/tabula-microcebus/kmermaid/dayhoff_ksize11/antoine__180917_A00111_0211_AHGKCVDMXX/ANTOINE_LUNG_P3/ANTOINE_LUNG_P3/ \
		--extract_coding_jaccard_threshold 0.5 \
		-work-dir ~/pureScratch/nextflow-intermediates/ \
		-resume
