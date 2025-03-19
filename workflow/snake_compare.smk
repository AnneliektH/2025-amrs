import os
import pandas as pd

# def ksizes
KSIZE = 21
# KSIZE = 25
# def outdir
OUTPUT_DIR ="/group/ctbrowngrp2/amhorst/2025-amrs/results/compare_all"
logs_dir = f"{OUTPUT_DIR}/logs"
ATLAS_DIR = "/group/ctbrowngrp2/scratch/annie/2023-swine-sra/results/atlas"
ALIGNERS = ['kma','bowtie2','bwa']

# set configfile
configfile: "../config/config.yaml" 
# Set samples for human and pig
metadata_metag = pd.read_csv(config['metag'], usecols=['acc'])
# Create a list of run ids
samples_pig = metadata_metag['acc'].tolist()
# Define samples
PIG_METAG = config.get('samples', samples_pig)


rule all:
    input:
        #expand("../results/compare_all/check/{metag}.rgi.{aligner}.done", metag=PIG_METAG, aligner=ALIGNERS),
        expand("../results/compare_all/{metag}.{ksize}.s1000.dbconcat.csv", metag=PIG_METAG, ksize=KSIZE),
        #expand("../results/compare_all/skipmer.{metag}.csv", metag=PIG_METAG)



# Use all
rule rgi:
    input:
        read_one=f"{ATLAS_DIR}/atlas_{{metag}}/{{metag}}/sequence_quality_control/{{metag}}_QC_R1.fastq.gz",
        read_two=f"{ATLAS_DIR}/atlas_{{metag}}/{{metag}}/sequence_quality_control/{{metag}}_QC_R2.fastq.gz",
    output:
        check = "../results/compare_all/check/{metag}.rgi.{aligner}.done",
    params:
        output_prefix=f"{OUTPUT_DIR}/rgi_{{metag}}.{{aligner}}"
    benchmark: f"{logs_dir}/rgi.{{metag}}.{{aligner}}.benchmark"
    conda: 
        "rgi"
    threads: 15
    shell: 
        """
        mkdir -p $(dirname {params.output_prefix})
        rgi bwt \
        -a {wildcards.aligner} \
        --read_one {input.read_one} \
        --read_two {input.read_two} \
        --output_file {params.output_prefix} \
        --local --clean -n {threads} && touch {output.check}
        """


# sourmash x the whole db
rule fastgather:
    input:
       sig = "/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/sketches_metag/{metag}.zip",
       db = '../resources/database/CARD_smash/card_protein_homolog.merge.nucleotide.zip'
    output:
        csv = "../results/compare_all/{metag}.{ksize}.s1000.dbconcat.csv",
    conda: 
        "branchwater-skipmer"
    threads: 15
    benchmark: f"{logs_dir}/fastgather.{{metag}}.{{ksize}}.benchmark"
    shell:
        """ 
        sourmash scripts fastgather \
        {input.sig} {input.db} \
        -k {wildcards.ksize} --scaled 1000 -t 0 \
        -c {threads} -o {output.csv}
        """

rule fastgather_skip:
    input:
       sig = "../resources/sketches/{metag}.skipmer.zip",
       db = '../resources/database/CARD_smash/nucleotide_card.prothomolog.skipmer.zip'
    output:
        csv = "../results/compare_all/skipmer.{metag}.csv",
    conda: 
        "branchwater-skipmer"
    threads: 15
    benchmark: f"{logs_dir}/fastgather.{{metag}}.skipmer.benchmark"
    shell:
        """ 
        sourmash scripts fastgather \
        {input.sig} {input.db} -m skipm2n3 -k 24 --scaled 100 -t 0 \
        -c {threads} -o {output.csv}
        """

