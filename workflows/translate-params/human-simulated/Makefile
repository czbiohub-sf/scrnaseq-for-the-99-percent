# --- RNAseq params --- #
GENCODE_DIR=/home/olga/data_lg/czbiohub-reference/gencode/human/v30
GTF=${GENCODE_DIR}/gencode.v30.annotation.ERCC92.gtf
FASTA=${GENCODE_DIR}/GRCh38.p12.genome.fa
STAR_INDEX=${GENCODE_DIR}/STARIndex


# --- Quest for orthologs data --- #
QFO_2019=/home/olga/data_sm/immune-evolution/databases/quest-for-orthologs/2019



# --- predictorthologs params -- #
PROTEOME_FASTA=${QFO_2019}/2019/Eukaryota/UP000005640_9606.fasta

# Fasta of the isoforms of the genes
ADDITIONAL_FASTA=${QFO_2019}/Eukaryota/UP000005640_9606_additional.fasta
COMBINED_FASTA=${QFO_2019}/Eukaryota/UP000005640_9606__combined_with_additional.fasta

MOUSE_FASTA=${QFO_2019}/Eukaryota/UP000000589_10090.fasta
BUSCO_MAMMALIA=/home/olga/data_sm/immune-evolution/databases/busco/orthodb-v10/mammalia_odb10/mammalia_odb10.fasta
BUSCO_MAMMALIA_HUMAN=${QFO_2019}/Eukaryota/UP000005640_9606__busco_mammlia_odbv10.fasta
BUSCO_MAMMALIA_MOUSE=${QFO_2019}/Eukaryota/UP000000589_10090__busco_mammlia_odbv10.fasta

OUTDIR_BASE=/mnt/data_sm/home/olga/pipeline-results/human-simulated/nf-predictorthologs
WORK_DIR=/mnt/data_sm/home/olga/nextflow-work/
READS_DIR=/home/olga/data_sm/immune-evolution/rawdata/simulated/human
READS_BUSCO_DIR=/home/olga/data_sm/immune-evolution/rawdata/simulated/human/busco_mammalia__orthodb-v10

COMMON_FLAGS=nextflow run czbiohub/nf-predictorthologs \
		-r olgabot/fix-skip-trimming \
		-latest \
		-profile docker \
		--skip_trimming \
		--skip_fastqc \
		--reads ${READS_DIR}/'*.fq.gz' \
		--proteome_search_fasta ${PROTEOME_FASTA} \
		--translate_jaccard_threshold 0.95 \
		--translate_peptide_ksize 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31 \
		--translate_peptide_molecule 'protein,dayhoff' \
		-with-tower \
		-w ${WORK_DIR} \
		--single_end \
		-resume \
		--max_cpus 90 \
		--max-memory 500.GB \
		--max_time 200.h \
		-dump-channels ch_translate_molecule_ksize \


BUSCO_COMMON_FLAGS=nextflow run czbiohub/nf-predictorthologs \
		-r olgabot/fix-skip-trimming \
		-latest \
		-profile docker \
		--skip_trimming \
		--skip_fastqc \
		--reads ${READS_BUSCO_DIR}/'*.fq.gz' \
		--proteome_search_fasta ${PROTEOME_FASTA} \
		--translate_jaccard_threshold 0.95 \
		--translate_peptide_ksize 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31 \
		--translate_peptide_molecule 'protein,dayhoff' \
		-with-tower \
		-w ${WORK_DIR} \
		--single_end \
		-resume \
		--max_cpus 90 \
		--max-memory 500.GB \
		--max_time 200.h


run:
	${COMMON_FLAGS} \
		--proteome_search_fasta ${PROTEOME_FASTA} \
		--proteome_translate_fasta ${PROTEOME_FASTA} \
		--outdir ${OUTDIR_BASE}--$@/


additional-isoforms:
	${COMMON_FLAGS} \
		--proteome_search_fasta ${COMBINED_FASTA} \
		--proteome_translate_fasta ${ADDITIONAL_FASTA} \
		--outdir ${OUTDIR_BASE}--$@/


run--search-additional:
	${COMMON_FLAGS} \
		--proteome_search_fasta ${COMBINED_FASTA} \
		--proteome_translate_fasta ${PROTEOME_FASTA} \
		--outdir ${OUTDIR_BASE}--$@/


busco-mammalia:
	${BUSCO_COMMON_FLAGS} \
		--proteome_search_fasta ${BUSCO_MAMMALIA} \
		--proteome_translate_fasta ${BUSCO_MAMMALIA} \
		--outdir ${OUTDIR_BASE}--$@/


busco-mammalia-human:
	${BUSCO_COMMON_FLAGS} \
		--proteome_search_fasta ${BUSCO_MAMMALIA} \
		--proteome_translate_fasta ${BUSCO_MAMMALIA_HUMAN} \
		--outdir ${OUTDIR_BASE}--$@/


busco-mammalia-mouse:
	${BUSCO_COMMON_FLAGS} \
		--proteome_search_fasta ${BUSCO_MAMMALIA} \
		--proteome_translate_fasta ${BUSCO_MAMMALIA_MOUSE} \
		--outdir ${OUTDIR_BASE}--$@/



rnaseq:
	nextflow run nf-core/rnaseq \
		-latest \
		-profile docker \
		-r 1.4.2 \
		--fasta ${FASTA} \
		--gtf ${GTF} \
		--star_index ${STAR_INDEX} \
		--skipTrimming \
		--skipFastQC \
		--reads ${READS_DIR}/'*.fq.gz' \
		-with-tower \
		-w ${WORK_DIR} \
		--singleEnd \
		-resume \
		--search_noncoding \
		--max_cpus 90 \
		--max-memory 500.GB \
		--max_time 200.h \
		--outdir ${OUTDIR_BASE}--$@/ \
		--fc_group_features_type gene_type

