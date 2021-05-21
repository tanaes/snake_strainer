from os.path import isfile, basename, dirname
from os import rename, remove
import shutil

from validator_collection.checkers import is_url
from urllib.request import urlretrieve

import gzip


def is_gzip(fp):
    with gzip.open(fp, 'r') as f:
        try:
            f.read(1)
            return(True)
        except OSError:
            return(False)


rule import_genome_list:
    """
    Imports list of genome bins for use.

    Takes an input file that is a list of paths and either
    copies (if local) or downloads (if URL) to 
    """
    input:
        genome_list=lambda wildcards: config['references'][wildcards.ref_set]['genome_list']
    output:
        genome_done='output/references/{ref_set}.txt'
    params:
        ref_prefix=lambda wildcards: config['references'][wildcards.ref_set]['prefix']
    threads: 1
    run:
        prefix = params.ref_prefix
        genome_dir = dirname(output.genome_done)
        path_list = ''
        print(input[0])
        genomes = pd.read_csv(input.genome_list,
                              sep='\t',
                              header=None,
                              index_col=0,
                              squeeze=True)

        for name, genome_path in genomes.iteritems():
            fname = join(genome_dir,
                         '{0}-{1}.fasta'.format(prefix,
                                                name))
            if isfile(genome_path):
                shutil.copy(genome_path, fname)
            elif is_url(genome_path):
                urlretrieve(genome_path, fname)

            if is_gzip(fname):
                rename(fname, '%s.gz' % fname)
                with gzip.open('%s.gz' % fname, 'rb') as f_in:
                    with open(fname, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove('%s.gz' % fname)

            path_list += '%s\n' % fname

        with open(output.genome_done, 'w') as f:
            f.write(path_list)


rule import_genomes:
    """
    Imports each genome list file
    """
    input:
        expand('output/references/{ref_set}.txt',
               ref_set=config['references'])
    output:
        'output/references.txt'
    run:
        with open(output[0], 'w') as o:
            for fp in input:
                with open(fp, 'r') as f:
                    for line in f:
                        o.write(line)

rule drep:
    """
    Runs dRep
    """
    input:
        'output/references.txt'
    output:
        drep_dir=directory('output/drep'),
        derep_genomes=directory('output/drep/dereplicated_genomes'),
        genome_info='output/drep/data_tables/genomeInformation.csv',
        prodigal=directory('output/drep/data/prodigal')
    conda:
        '../Envs/instrain.yaml'
    params:
        sa=config['params']['drep']['sa'],
        other=config['params']['drep']['other']
    log:
        'output/logs/drep.log'
    threads:
        config['params']['drep']['threads']
    benchmark:
        'output/benchmarks/drep.txt'
    shell:
        """
        dRep dereplicate {output.drep_dir} \
         -p {threads} \
         -g {input[0]} \
         -sa {params.sa} \
         {params.other} 2> {log} 1>&2
        """

