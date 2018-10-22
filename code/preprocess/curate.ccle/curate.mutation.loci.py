##################################################################################
###                         curate.mutation.loci.py                            ###
##################################################################################
# curate ccle mutation in loci level

proj_path = 'D:/projects/DLCell'
import os
import re
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util

geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
mutation_path = os.path.join(proj_path, 'data/DepMap/CCLE_Mutation_20180718.txt')
out_path = os.path.join(proj_path, 'data/curated/lung/ccle/lung_MutExpr_cancergene_loci.csv')


##############################      main      ################################
# read cancer gene list
genes = pd.read_csv(geneanno_path)['Gene'].tolist()

# read mutation data
mutation_all = pd.read_csv(mutation_path, sep='\t', usecols=['Hugo_Symbol', 'Variant_Classification', 'Genome_Change', 'Tumor_Sample_Barcode'])
mutation_all.columns = ['Gene', 'MutClass', 'Mut', 'Cell']

# filter for lung
mutations = mutation_all.loc[mutation_all.Cell.str.contains('LUNG'),:]
mutations['Cell'] = mutations.Cell.apply(lambda x: re.sub(r'_LUNG.*', '', x).upper())

# filter for cancer genes
mutations = mutations.loc[mutations.Gene.isin(genes),:]

# mutation info
mutations['MutChr'] = mutations['Mut'].apply(lambda x: re.findall(r'(chr.*):', x)[0])
mutations['MutStart'] = mutations['Mut'].apply(lambda x: re.findall(r':([0-9]*)', x)[0])
mutations['Mut'] = mutations['Mut'].apply(lambda x: re.findall(r':[0-9]*(.*)', x.replace('_', ''))[0])
mutations = mutations[['Cell', 'Gene', 'MutChr', 'MutStart', 'Mut', 'MutClass']]
mutations.sort_values(by=['Cell', 'MutStart', 'Gene'])
mutations.reset_index(drop=True, inplace=True)

# divide genes
gene_loc_encoder = util.geneLocEncoder(genes=mutations['Gene'], locs=mutations['MutStart'])
mutations_loc = util.divideGene(mutations, gene_loc_encoder)
mutations_loc.to_csv(out_path, index=None)
