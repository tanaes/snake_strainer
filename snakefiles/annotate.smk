import pandas as pd
from os.path import join
from os import makedirs
from collections import OrderedDict

localrules: rename_fasta, prep_bakta, annotate


def simplify_fasta(old_fp, new_fp, new_name):
    contig_name_dict = OrderedDict({})

    n = 0
    with open(old_fp, 'r') as f:
        with open(new_fp, 'w') as o:
            for line in f:
                if line.startswith('>'):
                    n += 1
                    contig_name = '{0}_{1}'.format(new_name,
                                                   n)
                    newline = '>{0} {1}\n'.format(contig_name,
                                                  line[1:].rstrip())
                    contig_name_dict['old_contig'] = contig_name
                else:
                    newline = line
                o.write(newline)

    return(contig_name_dict)


rule rename_fasta:
    input:
        genome_fna=lambda wildcards: genome_fps.loc[rev_rename_dict[wildcards.renamed], 'fp']
    output:
        renamed_fna='output/genomes/renamed/fastas/{renamed}.fna',
        contig_dict='output/genomes/renamed/contigs/{renamed}.txt'
    run:
        contig_names = simplify_fasta(input.genome_fna,
                                      output.renamed_fna,
                                      wildcards.renamed)
        with open(output.contig_dict, 'w') as f:
            for old in contig_names:
                f.write('{0}\t{1}\n'.format(old, contig_names[old]))

rule prep_bakta:
    output:
        touch('output/annotate/bakta/prep.done')

rule bakta:
    """
    Runs bakta
    """
    input:
        renamed_fna='output/genomes/renamed/fastas/{renamed}.fna'
    output:
        outfile='output/annotate/bakta/{renamed}/{renamed}.tsv',
        outdir=directory('output/annotate/bakta/{renamed}')
    params:
        db_path=config['params']['bakta']['db_path'],
        tmp_dir=config['params']['bakta']['tmp_dir']
    log:
        'output/logs/annotate/bakta/bakta-{renamed}.log'
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
        'output/benchmarks/bakta/bakta-{renamed}.log'
    shell:
        """
        bakta --db {params.db_path} \
            --output {output.outdir} \
            --prefix {wildcards.renamed} \
            --keep-contig-headers \
            --threads {threads} \
            --tmp-dir {params.tmp_dir} \
            {input.renamed_fna} \
            2> {log} 1>&2
        """

rule annotate:
    input:
        expand('output/annotate/bakta/{renamed}/{renamed}.tsv',
               renamed=genome_fps['renamed'])
