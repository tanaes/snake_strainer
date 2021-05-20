rule calc_stb:
    input:
        derep_genomes='output/drep/dereplicated_genomes'
    output:
        stb_file='output/instrain/input/dereplicated_genomes.stb',
        fasta_cat='output/instrain/input/dereplicated_genomes.fasta'
    conda:
        '../Envs/instrain.yaml'
    shell:
        """
        parse_stb.py --reverse \
         -f {input.derep_genomes}/* \
         -o {output.stb_file}

        cat {input.derep_genomes}/* > {output.fasta_cat}
        """


rule prep_drep:
    input:
        prodigal='output/drep/data/prodigal'
    output:
        fna_cat='output/instrain/input/genes/dereplicated_genomes.genes.fna',
        faa_cat='output/instrain/input/genes/dereplicated_genomes.genes.faa'
    shell:
        """
        cat {input.prodigal}/*.faa > {output.faa_cat}
        cat {input.prodigal}/*.fna > {output.fna_cat}
        """

rule instrain:
    input:
        rules.prep_drep.output,
        rules.calc_stb.output