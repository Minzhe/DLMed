##################################################################################
###                         curate.mutation.loci.py                            ###
##################################################################################
# curate utsw mike white mutation in loci level

proj_path = '/work/bioinformatics/s418336/projects/DLMed/'
import os
import re
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.join(proj_path, 'code'))
import matplotlib.pyplot as plt
import seaborn as sns
from utility import utility as util

genelist = 'allgene'

################################    cancer gene    ###################################
if genelist == 'cancergene':
    geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
    mutation_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.mutation.txt')
    out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_cancergene.csv')
    out_loci_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_cancergene_loci.csv')

    mutation = pd.read_table(mutation_path, sep='\t', usecols=['official', 'Cell', 'Position', 'Allele_Freq'])
    mutation.columns = ['Gene', 'Cell', 'MutStart', 'AF']
    mutation = mutation[['Cell', 'Gene', 'MutStart', 'AF']]

    # clean cell line name
    mutation.Cell = mutation.Cell.apply(util.cleanCellLine)

    # filter for cancer gene
    genes = pd.read_csv(geneanno_path)['Gene'].tolist()
    mutation = mutation.loc[mutation.Gene.isin(genes),:]
    mutation.sort_values(by=['Cell', 'Gene'], inplace=True)
    # sns.distplot(mutation.AF)
    # plt.savefig(os.path.join(proj_path, 'data/curated/Lung/utsw.mw/mutation_af.png'))
    mutation.to_csv(out_path, index=None)

    # divide genes
    gene_loc_encoder = util.geneLocEncoder(genes=mutation['Gene'], locs=mutation['MutStart'])
    mutation_loc = util.divideGene(mutation, gene_loc_encoder)

    mutation_loc.to_csv(out_loci_path, index=None)

################################    cancer gene    ###################################
elif genelist == 'allgene':
    mutation_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.mutation.txt')
    out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_allgene.csv')
    out_loci_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_allgene_loci_af.csv')

    mutation = pd.read_table(mutation_path, sep='\t', usecols=['official', 'Cell', 'Position', 'Allele_Freq'])
    mutation.columns = ['Gene', 'Cell', 'MutStart', 'AF']
    mutation = mutation[['Cell', 'Gene', 'MutStart', 'AF']]

    # clean cell line name
    mutation.Cell = mutation.Cell.apply(util.cleanCellLine)

    # filter for cancer gene
    mutation.sort_values(by=['Cell', 'Gene'], inplace=True)
    # sns.distplot(mutation.AF)
    # plt.savefig(os.path.join(proj_path, 'data/curated/Lung/utsw.mw/mutation_af.png'))
    mutation.to_csv(out_path, index=None)

    # divide genes
    gene_loc_encoder = util.geneLocEncoder(genes=mutation['Gene'], locs=mutation['MutStart'])
    mutation_loc = util.divideGene(mutation, gene_loc_encoder)

    mutation_loc.to_csv(out_loci_path, index=None)
