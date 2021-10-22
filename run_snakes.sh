

snakemake -s simulate_translate_assess.snakefile --configfiles config/QfO_vertebrates.yml --profile farm --cluster-config config/cluster_config.yml --jobs 5 -n
