#################################################################################
###                           curate.mutation.py                              ###
#################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import re
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.join(proj_dir, 'code'))
from utility import utility as util

data_path = os.path.join(proj_dir, 'data/PDX/PDX.mutation.txt')
gene_path = os.path.join(proj_dir, 'data/COSMIC/Census_all.csv')
out_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.mutation.cancergenes.csv')
out_loci_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.mutation.cancergenes.loci.csv')

##############################    main   ################################
cancergenes = list(pd.read_csv(gene_path, usecols=[0], squeeze=True))
mut = pd.read_table(data_path, sep='\t')
mut = mut.loc[mut.Category.str.contains('Mut'),['Sample', 'Gene', 'Details']]
mut.Details = mut.Details.apply(lambda x: x.split(',')[0].replace('?', ''))
mut.Details = mut.Details.apply(lambda x: int(re.findall(r'-?[0-9]+', x)[0]))
mut.columns = ['Sample', 'Gene', 'Site']
mut = mut.loc[mut.Gene.isin(cancergenes),:]
mut = mut.loc[~mut.duplicated(),:]
mut.to_csv(out_path, index=None)

# divide genes
gene_loc_encoder = util.geneLocEncoder(genes=mut['Gene'], locs=mut['Site'], maxwin=96)
mut.columns = ['Cell', 'Gene', 'MutStart']
mut_loc = util.divideGene(mut, gene_loc_encoder)
mut_loc.columns = ['Sample', 'Gene', 'Loci']
mut_loc.to_csv(out_loci_path, index=None)
