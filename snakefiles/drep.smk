from os.path import isfile, basename, dirname
import shutil

from validator_collection.checkers import is_url
from urllib.request import urlretrieve

rule import_genome_list:
    """
    Imports list of genome bins for use.

    Takes an input file that is a list of paths and either
    copies (if local) or downloads (if URL) to 
    """
    input:
        genome_list=lambda wildcards: config['references'][wildcards.ref_set]['genome_list']
    output:
        genome_done=touch('output/references/{ref_set}.done')
    params:
        ref_prefix=lambda wildcards: config['references'][wildcards.ref_set]['prefix']
    threads: 1
    run:
        prefix = params.ref_prefix
        genome_dir = dirname(output.genome_done)
        print(input[0])
        genomes = pd.read_csv(input.genome_list,
                              sep='\t',
                              header=None,
                              index_col=0,
                              squeeze=True)

        for name, genome_path in genomes.iteritems():
            fname = join(genome_dir,
                         '{0}-{1}.fastq.gz'.format(prefix,
                                                   name))
            if isfile(genome_path):
                shutil.copy(genome_path, fname)
            elif is_url(genome_path):
                urlretrieve(genome_path, fname)

rule import_genomes:
    """
    Imports each genome list file
    """
    input:
        expand('output/references/{ref_set}.done',
               ref_set=config['references'])
    output:
        touch('output/references.done')

