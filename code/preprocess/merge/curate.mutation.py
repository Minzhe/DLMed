##################################################################################
###                         curate.mutation.loci.py                            ###
##################################################################################
# curate utsw mike white mutation in loci level

proj_path = 'D:/projects/DLMed'
import os
import re
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util

mutation_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_Mutation_cancergene.csv')
out_loci_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_Mutation_cancergene_loci.csv')

##############################      main      ################################
mutation = pd.read_csv(mutation_path)

# divide genes
gene_loc_encoder = util.geneLocEncoder(genes=mutation['Gene'], locs=mutation['MutStart'])
mutation_loc = util.divideGene(mutation, gene_loc_encoder)
mutation_loc.to_csv(out_loci_path, index=False)