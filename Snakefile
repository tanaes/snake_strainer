import pandas as pd
from os.path import join, splitext
from yaml import safe_load
from os import makedirs

def get_read(sample, read):
    return(samples_df.loc[sample, read])

def input_genomes(ref_fp, prefix='sGB'):
    genomes = pd.read_csv(ref_fp, 
                          sep='\t', 
                          header=0, 
                          index_col=0)

    if 'fp' not in genomes.columns:
        raise ValueError('must have fp column to fasta file')

    if 'renamed' not in genomes.columns:
        new_names = ['{0}_{1}'.format(prefix,
                                      x)
                     for x in range(genomes.shape[0])]

        genomes['renamed'] = new_names

        genomes.to_csv(splitext(ref_fp)[0] + '.renamed.' + splitext(ref_fp)[1])

    return(genomes)


configfile: 'config.yaml'

with open(config['resources'], 'r') as f:
    res = safe_load(f)

references = config['references']
samples_df = pd.read_csv(config['samples_fp'],
                         sep='\t',
                         header=0,
                         index_col=0)

samples = list(samples_df.index)
prefix=config['params']['simplify_fastas']['prefix']

input_fps = input_genomes(references['genome_list'])

genome_fps = input_fps['fp']

print(genome_fps)
genomes = list(genome_fps.index)
print(genomes)

rev_rename_dict = {}
for x in genome_fps.index:
    rev_rename_dict[genome_fps.loc[x, 'renamed']] = x

include: 'snakefiles/instrain.smk'
include: 'snakefiles/annotate.smk'


rule all:
    input:
        rules.prep_drep.output,
        rules.instrain_compare.output
