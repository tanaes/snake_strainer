localrules: prep_drep_input, prep_drep


rule prep_drep_input:
    input:
        expand('output/genomes/renamed/fastas/{renamed}.fna',
               renamed=genome_fps['renamed'])
    output:
        'output/instrain/input/bakta_genome_list.txt'
    run:
        with open(output[0], 'w') as f:
            for fp in input:
                f.write('%s\n' % fp)

# rule prep_drep_prodigal:
#     output:
#         fna_cat='output/instrain/input/prodigal/genes/dereplicated_genomes.genes.fna',
#         faa_cat='output/instrain/input/prodigal/genes/dereplicated_genomes.genes.faa',
#         stb_file='output/instrain/input/prodigal/dereplicated_genomes.stb',
#         fasta_cat='output/instrain/input/prodigal/dereplicated_genomes.fasta'
#     params:
#         derep_genomes=config['references']['drep']['prodigal'],
#         prodigal=config['references']['drep']['derep_genomes'],
#     conda:
#         '../Envs/instrain.yaml'
#     shell:
#         """
#         cat {params.prodigal}/*.faa > {output.faa_cat}
#         cat {params.prodigal}/*.fna > {output.fna_cat}

#         parse_stb.py --reverse \
#          -f {params.derep_genomes}/* \
#          -o {output.stb_file}

#         cat {params.derep_genomes}/* > {output.fasta_cat}
#         """

rule cat_fasta:
    input:
        fastas=lambda wildcards: expand('output/annotate/bakta/{renamed}/{renamed}.{ext}',
                                        renamed=genome_fps['renamed'],
                                        ext=wildcards.ext)
    output:
        fasta_fp='output/instrain/input/dereplicated_genomes.{ext}'
    run:
        with open(output[0], 'w') as o:
            for fp in input.fastas:
                with open(fp, 'r') as f:
                    for line in f:
                        o.write(line)

rule cat_fastas:
    input:
        expand('output/instrain/input/dereplicated_genomes.{ext}',
               ext=['fna','faa','ffn','gbff'])

rule prep_stb_bakta:
    input:
        fasta_fps='output/instrain/input/bakta_genome_list.txt'
    output:
        stb_file='output/instrain/input/bakta/dereplicated_genomes.stb'
    params:
        derep_genomes=config['references']['drep']['prodigal'],
        prodigal=config['references']['drep']['derep_genomes'],
    conda:
        '../Envs/instrain.yaml'
    shell:
        """
        parse_stb.py --reverse \
         -f {input.fasta_fps} \
         -o {output.stb_file}
        """

rule index_db:
    input:
        reference='output/instrain/input/dereplicated_genomes.fna'
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
    benchmark:
        "output/benchmarks/instrain/index_db"
    shell:
        """
        bowtie2-build --threads {threads} {params.other} \
        {input.reference} output/instrain/input/references 2> {log} 1>&2
        """

rule map_reads:
    input:
        fwd=join(config['reads_dir'], '{sample}.R1.fastq.gz'),
        rev=join(config['reads_dir'], '{sample}.R2.fastq.gz'),
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
    benchmark:
        "output/benchmarks/instrain/map_reads/map_reads.{sample}"
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

def get_mem_mb(wildcards, attempt, default=4000):
    return attempt * default

rule instrain_profile:
    input:
        aln=rules.map_reads.output.aln,
        reference='output/instrain/input/dereplicated_genomes.fna',
        genes_file='output/instrain/input/dereplicated_genomes.gbff',
        stb_file='output/instrain/input/bakta/dereplicated_genomes.stb'
    output:
        profile=directory('output/instrain/output/profiles/{sample}.IS'),
        bam='output/instrain/input/alignments/{sample}.sorted.bam'
    conda:
        "../Envs/instrain.yaml"
    log:
        "output/logs/instrain/instrain_profile.{sample}.log"
    threads:
        res['instrain_profile']['threads']
    resources:
        partition = res['instrain_profile']['partition'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(wildcards, attempt, default=res['instrain_profile']['mem_mb']),
        qos = res['instrain_profile']['qos'],
        time = res['instrain_profile']['time']
    benchmark:
        "output/benchmarks/instrain/profile/profile.{sample}"
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
        --database_mode  2> {log} 1>&2
        """


rule instrain_compare:
    input:
        profiles=expand(rules.instrain_profile.output.profile,
                        sample=samples),
        stb_file='output/instrain/input/bakta/dereplicated_genomes.stb'
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
    benchmark:
        "output/benchmarks/instrain/compare"
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
