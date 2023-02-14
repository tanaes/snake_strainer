localrules: prep_drep

rule prep_drep:
    output:
        fna_cat='output/instrain/input/genes/dereplicated_genomes.genes.fna',
        faa_cat='output/instrain/input/genes/dereplicated_genomes.genes.faa',
        stb_file='output/instrain/input/dereplicated_genomes.stb',
        fasta_cat='output/instrain/input/dereplicated_genomes.fasta'
    params:
        derep_genomes=config['references']['drep']['prodigal'],
        prodigal=config['references']['drep']['derep_genomes'],
    conda:
        '../Envs/instrain.yaml'
    shell:
        """
        cat {params.prodigal}/*.faa > {output.faa_cat}
        cat {params.prodigal}/*.fna > {output.fna_cat}

        parse_stb.py --reverse \
         -f {params.derep_genomes}/* \
         -o {output.stb_file}

        cat {params.derep_genomes}/* > {output.fasta_cat}
        """


rule index_db:
    input:
        reference=rules.prep_drep.output.fasta_cat
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
        res['bowtie2_build']['threads']
    resources:
        partition = res['bowtie2_build']['partition'],
        mem_mb = res['bowtie2_build']['mem_mb'],
        qos = res['bowtie2_build']['qos'],
        time = res['bowtie2_build']['time']
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
        res['map_reads']['threads']
    resources:
        partition = res['map_reads']['partition'],
        mem_mb = res['map_reads']['mem_mb'],
        qos = res['map_reads']['qos'],
        time = res['map_reads']['time']
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
        reference=rules.prep_drep.output.fasta_cat,
        genes_file=rules.prep_drep.output.fna_cat,
        stb_file=rules.prep_drep.output.stb_file
    output:
        profile=directory('output/instrain/output/profiles/{sample}.IS'),
        bam='output/instrain/input/alignments/{sample}.sorted.bam'
    conda:
        "../Envs/instrain.yaml"
    threads:
        res['instrain_profile']['threads']
    resources:
        partition = res['instrain_profile']['partition'],
        mem_mb = res['instrain_profile']['mem_mb'],
        qos = res['instrain_profile']['qos'],
        time = res['instrain_profile']['time']
    params:
        other=config['params']['instrain']['profile']
    shell:
        """
        inStrain profile {input.aln} \
        {input.reference} \
        -o {output.profile} \
        -p {threads} \
        -g {input.genes_file} \
        -s {input.stb_file} \
        {params.other} \
        --database_mode
        """


rule instrain_compare:
    input:
        profiles=expand(rules.instrain_profile.output.profile,
                        sample=samples),
        stb_file=rules.prep_drep.output.stb_file
    output:
        compare=directory('output/instrain/output/compare')
    threads:
        res['instrain_compare']['threads']
    resources:
        partition = res['instrain_compare']['partition'],
        mem_mb = res['instrain_compare']['mem_mb'],
        qos = res['instrain_compare']['qos'],
        time = res['instrain_compare']['time']
    params:
        other=config['params']['instrain']['compare']
    conda:
        "../Envs/instrain.yaml"
    log:
        "output/logs/instrain/instrain_compare.log"
    shell:
        """
        inStrain compare \
        -s {input.stb_file} \
        -p {threads} \
        --database_mode \
        -o {output.compare} \
        -i {input.profiles} \
        {params.other} \
        2> {log} 1>&2

        """


rule coverage_calc:
    input:
        bam=rules.instrain_profile.output.bam
    output:
        cov='output/instrain/input/alignments/{sample}.cov'
    conda:
        "../Envs/instrain.yaml"
    log:
        "output/logs/instrain/{sample}.coverage_calc.log"
    shell:
        """
        bedtools genomecov -ibam {input.bam} -dz > {output.cov} 2> {log}
        """


rule coverage:
    input:
        expand(rules.coverage_calc.output.cov,
               sample=samples)
