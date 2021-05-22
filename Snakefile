import pandas as pd
from os.path import join

configfile: 'config.yaml'

include: 'snakefiles/drep.smk'
include: 'snakefiles/instrain.smk'

references = config['references']
samples_df = pd.read_csv(config['samples_fp'],
                         sep='\t',
                         header=0,
                         index_col=0)

samples = list(samples_df.index)


def get_read(sample, read):
    return(samples_df.loc[sample, read])


rule all:
    input:
        rules.import_genomes.output,
        rules.drep.output,
        rules.prep_drep.output,
        rules.calc_stb.output,
        rules.instrain_compare.output
