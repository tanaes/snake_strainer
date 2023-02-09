import pandas as pd
from os.path import join
from yaml import safe_load
from os import makedirs

def get_read(sample, read):
    return(samples_df.loc[sample, read])

def input_genomes(ref_fp):
    genomes = pd.read_csv(ref_fp, 
                          sep='\t', 
                          header=None, 
                          index_col=0)
    return(genomes[1])

def simplify_fastas(genome_fps, prefix='sGB'):
    genome_table_fp = 'output/genomes/genome_translation.csv'
    contig_table_fp = 'output/genomes/contig_translation.csv'
    rename_dir = 'output/genomes/renamed'

    makedirs(rename_dir, exist_ok=True)
    genome_name_dict = {'old_genome': [],
                        'new_genome': []}
    contig_name_dict = {'old_genome': [],
                        'new_genome': [],
                        'old_contig': [],
                        'new_contig': []}

    renamed_fps = {}

    for i, genome in enumerate(genome_fps.index):
        genome_fp = genome_fps[genome]
        new_name = '{0}_{1:05d}'.format(prefix,
                                        i+1)
        genome_name_dict['old_genome'].append(genome)
        genome_name_dict['new_genome'].append(new_name)
        
        new_fp = join(rename_dir, new_name + '.fna')
        
        n = 0
        with open(genome_fp, 'r') as f:
            with open(new_fp, 'w') as o:
                for line in f:
                    if line.startswith('>'):
                        n += 1
                        contig_name = '{0}_{1}'.format(new_name,
                                                       n)
                        newline = '>{0} {1}\n'.format(contig_name,
                                                      line[1:].rstrip())
                        contig_name_dict['old_genome'].append(genome)
                        contig_name_dict['new_genome'].append(new_name)
                        contig_name_dict['old_contig'].append(line[1:].rstrip().split(' ')[0])
                        contig_name_dict['new_contig'].append(contig_name)
                    else:
                        newline = line
                    o.write(newline)
        
        pd.DataFrame(genome_name_dict).to_csv(genome_table_fp, 
                                              index=False)
        
        pd.DataFrame(contig_name_dict).to_csv(contig_table_fp, 
                                              index=False)

        renamed_fps[new_name] = new_fp

    return(pd.Series(renamed_fps))


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
genome_fps = simplify_fastas(input_fps, 
                             prefix)

print(genome_fps)
genomes = list(genome_fps.index)
print(genomes)

include: 'snakefiles/instrain.smk'
include: 'snakefiles/annotate.smk'


rule all:
    input:
        rules.prep_drep.output,
        rules.instrain_compare.output
