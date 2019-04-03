##################################################################################
###                         curate.mutation.loci.py                            ###
##################################################################################
# curate ccle mutation in loci level

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

genelist = 'allgenes'

##############################      cancer genes      ################################
if genelist == 'cancergenes':
    geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
    mutation_path = os.path.join(proj_path, 'data/DepMap/CCLE_Mutation_20180718.txt')
    out_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_Mutation_cancergene.csv')
    out_loci_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_Mutation_cancergene_loci_af.csv')

    # read cancer gene list
    genes = pd.read_csv(geneanno_path)['Gene'].tolist()

    # read mutation data
    mutation_all = pd.read_csv(mutation_path, sep='\t', usecols=['Hugo_Symbol', 'Variant_Classification', 'Genome_Change', 'Tumor_Sample_Barcode', 'WES_AC'])
    mutation_all.columns = ['Gene', 'MutClass', 'Mut', 'Cell', 'AF']

    # allele frequency
    mutation_all.AF = mutation_all.AF.apply(str)
    af0 = mutation_all.AF.apply(lambda x: x.split(':')[0]).astype('float')
    af1 = mutation_all.AF.apply(lambda x: x.split(':')[-1]).astype('float')
    mutation_all.AF = af0 / (af0 + af1)
    sns.distplot(mutation_all.AF.dropna())
    plt.savefig(os.path.join(proj_path, 'data/curated/Lung/ccle/mutation_af.png'))

    # filter for lung
    mutations = mutation_all.loc[mutation_all.Cell.str.contains('LUNG'),:]
    mutations['Cell'] = mutations.Cell.apply(lambda x: re.sub(r'_LUNG.*', '', x).upper())

    # filter for cancer genes
    mutations = mutations.loc[mutations.Gene.isin(genes),:]

    # mutation info
    mutations['MutChr'] = mutations['Mut'].apply(lambda x: re.findall(r'(chr.*):', x)[0])
    mutations['MutStart'] = mutations['Mut'].apply(lambda x: re.findall(r':([0-9]*)', x)[0])
    mutations['Mut'] = mutations['Mut'].apply(lambda x: re.findall(r':[0-9]*(.*)', x.replace('_', ''))[0])
    mutations = mutations[['Cell', 'Gene', 'MutChr', 'MutStart', 'Mut', 'MutClass', 'AF']]
    mutations.AF = mutations.AF.apply(lambda x: 0.5 if np.isnan(x) else x)
    mutations.sort_values(by=['Cell', 'MutStart', 'Gene'])
    mutations.reset_index(drop=True, inplace=True)
    mutations.to_csv(out_path, index=None)

    # divide genes
    gene_loc_encoder = util.geneLocEncoder(genes=mutations['Gene'], locs=mutations['MutStart'])
    mutations_loc = util.divideGene(mutations, gene_loc_encoder)
    mutations_loc.to_csv(out_loci_path, index=None)



##############################      all genes      ################################
elif genelist == 'allgenes':
    mutation_path = os.path.join(proj_path, 'data/DepMap/CCLE_Mutation_20180718.txt')
    out_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_Mutation_allgene.csv')
    out_loci_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_Mutation_allgene_loci_af.csv')

    # read mutation data
    mutation_all = pd.read_csv(mutation_path, sep='\t', usecols=['Hugo_Symbol', 'Variant_Classification', 'Genome_Change', 'Tumor_Sample_Barcode', 'WES_AC'])
    mutation_all.columns = ['Gene', 'MutClass', 'Mut', 'Cell', 'AF']

    # allele frequency
    mutation_all.AF = mutation_all.AF.apply(str)
    af0 = mutation_all.AF.apply(lambda x: x.split(':')[0]).astype('float')
    af1 = mutation_all.AF.apply(lambda x: x.split(':')[-1]).astype('float')
    mutation_all.AF = af0 / (af0 + af1)
    sns.distplot(mutation_all.AF.dropna())
    plt.savefig(os.path.join(proj_path, 'data/curated/Lung/ccle/mutation_allgene_af.png'))

    # filter for lung
    mutations = mutation_all.loc[mutation_all.Cell.str.contains('LUNG'),:]
    mutations['Cell'] = mutations.Cell.apply(lambda x: re.sub(r'_LUNG.*', '', x).upper())

    # mutation info
    mutations['MutChr'] = mutations['Mut'].apply(lambda x: re.findall(r'(chr.*):', x)[0])
    mutations['MutStart'] = mutations['Mut'].apply(lambda x: re.findall(r':([0-9]*)', x)[0])
    mutations['Mut'] = mutations['Mut'].apply(lambda x: re.findall(r':[0-9]*(.*)', x.replace('_', ''))[0])
    mutations = mutations[['Cell', 'Gene', 'MutChr', 'MutStart', 'Mut', 'MutClass', 'AF']]
    mutations = mutations.loc[mutations.MutClass != 'Silent',:]
    mutations.AF = mutations.AF.apply(lambda x: 0.5 if np.isnan(x) else x)
    mutations.sort_values(by=['Cell', 'MutStart', 'Gene'])
    mutations.reset_index(drop=True, inplace=True)
    mutations.to_csv(out_path, index=None)

    # divide genes
    gene_loc_encoder = util.geneLocEncoder(genes=mutations['Gene'], locs=mutations['MutStart'])
    mutations_loc = util.divideGene(mutations, gene_loc_encoder)
    mutations_loc.to_csv(out_loci_path, index=None)
