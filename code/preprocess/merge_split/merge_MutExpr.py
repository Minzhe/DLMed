#################################################################################
###                             merge_MutExpr.py                              ###
#################################################################################
# merge mutation allele frequency and expression data

proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os, sys
import numpy as np
import pandas as pd
sys.path.append(os.path.join(proj_dir, 'code'))
import matplotlib.pyplot as plt
plt.style.use('seaborn')
import seaborn as sns
import utility.plot as p
import utility.integrate as itg
import utility.utility as util


#########################        main        #########################
geneanno_path = os.path.join(proj_dir, 'data/curated/cancer.gene.anno.csv')
# mutation
ccle_mut_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/ccle.lung_Mutation_cancergene.csv')
utsw_mut_path = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_cancergene.csv')
mut_out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mutation_cancergene.csv')
out_loci_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mutation_cancergene.loci.csv')

ccle_mutation = pd.read_csv(ccle_mut_path)
utsw_mutation = pd.read_csv(utsw_mut_path)
genes = pd.read_csv(geneanno_path)['Gene'].tolist()
ccle_mutation = ccle_mutation.loc[ccle_mutation.Gene.isin(genes),:]
utsw_mutation = utsw_mutation.loc[utsw_mutation.Gene.isin(genes),:]

mutation = itg.mergeMutation(CCLE=ccle_mutation, UTSW=utsw_mutation)
mutation.to_csv(mut_out_path, index=None)
# sns.distplot(mutation.AF[mutation.AF!=0.5].dropna())
# plt.savefig(os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/mutation_af.png'))

# divide genes
gene_loc_encoder = util.geneLocEncoder(genes=mutation['Gene'], locs=mutation['MutStart'])
mutation_loc = util.divideGene(mutation, gene_loc_encoder)
mutation_loc.to_csv(out_loci_path, index=False)

# expression
ccle_expr_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/ccle.lung_RNAseq_cancergene.csv')
utsw_expr_path = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/utsw.lung_RNAseq_cancergene.csv')
expr_out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_RNAseq_cancergene.csv')

ccle_expr = pd.read_csv(ccle_expr_path, index_col=0)
utsw_expr = pd.read_csv(utsw_expr_path, index_col=0)

expr = itg.mergeExpression(CCLE=ccle_expr, UTSW=utsw_expr, keep_dup=['CCLE', 'UTSW'])
expr.to_csv(expr_out_path)