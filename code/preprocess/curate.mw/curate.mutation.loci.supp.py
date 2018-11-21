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

geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
mutation_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.mutation.txt')
out_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_Mutation_cancergene.csv')
out_loci_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_Mutation_cancergene_loci.csv')

################################    main    ###################################
mutation = pd.read_table(mutation_path, sep='\t', usecols=['official', 'Cell', 'Position'])
mutation.columns = ['Gene', 'Cell', 'MutStart']
mutation = mutation[['Cell', 'Gene', 'MutStart']]

# clean cell line name
mutation.Cell = mutation.Cell.apply(util.cleanCellLine)

# filter for cancer gene
genes = pd.read_csv(geneanno_path)['Gene'].tolist()
mutation = mutation.loc[mutation.Gene.isin(genes),:]
mutation.sort_values(by=['Cell', 'Gene'], inplace=True)
mutation.to_csv(out_path, index=None)

# divide genes
gene_loc_encoder = util.geneLocEncoder(genes=mutation['Gene'], locs=mutation['MutStart'])
mutation_loc = util.divideGene(mutation, gene_loc_encoder)

mutation_loc.to_csv(out_loci_path, index=None)