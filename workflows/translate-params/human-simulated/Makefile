PROTEOME_FASTA=/home/olga/data_sm/kmer-hashing/quest-for-orthologs/data/2019/Eukaryota/UP000005640_9606.fasta
OUTDIR_BASE=/mnt/data_sm/home/olga/pipeline-results/human-simulated/nf-predictorthologs
WORK_DIR=/mnt/data_sm/home/olga/nextflow-work/
READS_DIR=/home/olga/data_sm/immune-evolution/rawdata/simulated/human


run:
	nextflow run czbiohub/nf-predictorthologs \
		-r  olgabot/fix-skip-trimming \
		-profile docker \
		--skip_trimming \
		--skip_fastqc \
		--reads ${READS_DIR}/'*.fq.gz' \
		--proteome_search_fasta ${PROTEOME_FASTA} \
		--proteome_translate_fasta ${PROTEOME_FASTA} \
		--translate_jaccard_threshold 0.95 \
		--translate_peptide_ksize 9,10,11,12,13,14,15,16,17,18,19,20,21 \
		--translate_molecule protein,dayhoff \
		-with-tower \
		--outdir ${OUTDIR_BASE}--$@/ \
		-w ${WORK_DIR} \
		--single_end \
		-resume \
		--search_noncoding \
		--max_cpus 90 \
		--max-memory 500.GB \
		--max_time 200.h 
