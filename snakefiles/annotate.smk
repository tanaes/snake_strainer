import pandas as pd
from os.path import join


rule prep_bakta:
    output:
        touch('output/annotate/bakta/prep.done')


rule bakta:
    """
    Runs bakta
    """
    input:
        genome_fna=lambda wildcards: genome_fps[wildcards.genome]
    output:
        outfile='output/annotate/bakta/{genome}/{genome}.tsv',
        outdir='output/annotate/bakta/{genome}'
    params:
        db_path=config['params']['bakta']['db_path'],
        tmp_dir=config['params']['bakta']['tmp_dir']
    log:
        'output/logs/annotate/bakta/bakta-{genome}.log'
    threads:
        res['bakta']['threads']
    conda:
        '../Envs/bakta.yaml'
    resources:
        partition = res['bakta']['partition'],
        mem_mb = res['bakta']['mem_mb'],
        qos = res['bakta']['qos'],
        time = res['bakta']['time']
    benchmark:
        'output/benchmarks/bakta/bakta-{genome}.log'
    shell:
        """
        bakta --db {params.db_path} \
            --verbose \
            --output {output.outdir} \
            --prefix {wildcards.genome} \
            --threads {threads} \
            --tmp-dir {params.tmp_dir}
            {input.genome_fna} \
            2> {log} 1>&2
        """

rule annotate:
    input:
        expand('output/annotate/bakta/{genome}/{genome}.tsv',
               genome=genomes)

