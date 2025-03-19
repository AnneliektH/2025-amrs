FASTA, = glob_wildcards('/group/ctbrowngrp2/scratch/annie/2024-pigparadigm/results/MAGs/genomes/{fasta}.fasta')


rule all:
    input:
        expand("../results/rgi_genomes/check/{genome}.done", genome=FASTA),
       

# Use all
rule rgi:
    input:
        contig='/group/ctbrowngrp2/scratch/annie/2024-pigparadigm/results/MAGs/genomes/{genome}.fasta'
    output:
        check = "../results/rgi_genomes/check/{genome}.done",
    conda: 
        "rgi"
    threads: 10
    shell: 
        """
        rgi main \
        -i {input.contig} -o ../results/rgi_genomes/{wildcards.genome} \
        -t contig -a DIAMOND -n 12 --clean \
        --keep && touch {output.check}
        """
