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


rule index_db:
    input:
        reference=rules.calc_stb.output.fasta_cat
    output:
        multiext('output/instrain/input/references',
                 ".1.bt2",
                 ".2.bt2",
                 ".3.bt2",
                 ".4.bt2",
                 ".rev.1.bt2",
                 ".rev.2.bt2")
    log:
        "output/logs/index_db.log"
    conda:
        "../Envs/bowtie2.yaml"
    threads:
        config['params']['bowtie2']['threads']
    params:
        other=config['params']['bowtie2']['index']
    shell:
        """
        bowtie2-build --threads {threads} {params.other} \
        {input.reference} output/instrain/input/references 2> {log} 1>&2
        """


rule map_reads:
    input:
        fwd=lambda wildcards: get_read(wildcards.sample,
                                       'R1'),
        rev=lambda wildcards: get_read(wildcards.sample,
                                       'R2'),
        db=rules.index_db.output
    output:
        aln='output/instrain/input/alignments/{sample}.sam'
    conda:
        "../Envs/bowtie2.yaml"
    log:
        "output/logs/map_reads/map_reads-{sample}.log"
    threads:
        config['params']['bowtie2']['threads']
    params:
        other=config['params']['bowtie2']['map']
    shell:
        """
        bowtie2 -p {threads} \
        -x output/instrain/input/references \
        -1 {input.fwd} \
        -2 {input.rev} \
        {params.other} \
        > {output.aln} 2> {log}
        """



rule instrain_profile:
    input:
        aln=rules.map_reads.output.aln,
        reference=rules.calc_stb.output.fasta_cat,
        fna_cat=rules.prep_drep.output.fna_cat,
        faa_cat=rules.prep_drep.output.faa_cat,
        stb_file=rules.calc_stb.output.stb_file
    output:
        profile='output/instrain/output/profiles/{sample}.IS'
    conda:
        "../Envs/instrain.yaml"
    threads: 4
    shell:
        """
        inStrain profile {input.aln} \
        {input.reference} \
        -o {output.profile} \
        -p {threads} \
        -g {input.fna_cat} \
        -s {input.stb_file} \
        --database_mode
        """
