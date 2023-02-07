import pandas as pd
from os.path import join
from yaml import safe_load


def get_read(sample, read):
    return(samples_df.loc[sample, read])

def input_genomes(ref_fp):
    genomes = pd.read_csv(ref_fp, 
                          sep='\t', 
                          header=None, 
                          index_col=0)
    return(genomes[1])

configfile: 'config.yaml'

with open(config['resources'], 'r') as f:
    res = safe_load(f)

references = config['references']
samples_df = pd.read_csv(config['samples_fp'],
                         sep='\t',
                         header=0,
                         index_col=0)

samples = list(samples_df.index)

genome_fps = input_genomes(references['genome_list'])
genomes = list(genome_fps.index)

include: 'snakefiles/instrain.smk'
include: 'snakefiles/annotate.smk'


rule all:
    input:
        rules.prep_drep.output,
        rules.instrain_compare.output
