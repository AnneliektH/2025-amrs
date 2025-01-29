
# subset to only things found by either sourmash or RGI
rule subset_fasta:
    input:
        subset="../results/compare_all/found_both.txt",
        db="/resources/database/CARD/nucleotide_fasta_protein_homolog_model.fasta",
    output:
        subfasta = "../results/compare_all/found_both.fa",
    conda: 
        "bbmap"
    threads: 1
    shell: 
        """
        filterbyname.sh in={input.db} out={output.subfasta} names={input.subset} \
        include=t substring=t
        """

rule kma_db:
    input:
        subset="../results/compare_all/found_both.fa",
    output:
        check = "../results/compare_all/check/found_both.kma.txt"
    params:
        output_folder="../results/compare_all/found_both.kma"
    conda: 
        "bowtie2"
    threads: 4
    shell: 
        """
        kma index -i {input.subset} -o {params.output_folder} \
        && touch {output.check}
        """

rule bowtie_db:
    input:
        subset="../results/compare_all/found_both.fa",
    output:
        check = "../results/compare_all/check/found_both.bt2.txt"
    params:
        output_prefix="../results/compare_all/found_both.bt2"
    conda: 
        "bowtie2"
    threads: 4
    shell: 
        """
        bowtie2-build {input.subset} {params.output_prefix} \
        -p {threads} && touch {output.check}
        """