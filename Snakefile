import pandas as pd
from os.path import join

configfile: "config.yaml"

include: "snakefiles/drep.smk"

references = config['references']

rule all:
	input:
		rules.import_genomes.output