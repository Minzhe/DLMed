#####################################################################################
###                             merge_MutLociExpr.py                              ###
#####################################################################################
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

genelist = 'allgene'

#########################        cancer gene        #########################
if genelist == 'cancergene':
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


#########################        all gene        #########################
elif genelist == 'allgene':
    # mutation
    ccle_mut_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/ccle.lung_Mutation_allgene.csv')
    utsw_mut_path = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_allgene.csv')
    mut_out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mutation_allgene.csv')
    out_loci_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mutation_allgene_loci_af.csv')

    ccle_mutation = pd.read_csv(ccle_mut_path)
    utsw_mutation = pd.read_csv(utsw_mut_path)

    mutation = itg.mergeMutation(CCLE=ccle_mutation, UTSW=utsw_mutation)
    mutation.to_csv(mut_out_path, index=None)
    # sns.distplot(mutation.AF[mutation.AF!=0.5].dropna())
    # plt.savefig(os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/mutation_af.png'))

    # divide genes
    gene_loc_encoder = util.geneLocEncoder(genes=mutation['Gene'], locs=mutation['MutStart'])
    mutation_loc = util.divideGene(mutation, gene_loc_encoder)
    mutation_loc.to_csv(out_loci_path, index=False)

    # expression
    ccle_expr_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/ccle.lung_RNAseq_allgene.csv')
    utsw_expr_path = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/utsw.lung_RNAseq_allgene.csv')
    expr_out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_RNAseq_allgene.csv')

    ccle_expr = pd.read_csv(ccle_expr_path, index_col=0)
    utsw_expr = pd.read_csv(utsw_expr_path, index_col=0)

    expr = itg.mergeExpression(CCLE=ccle_expr, UTSW=utsw_expr, keep_dup=['CCLE', 'UTSW'])
    expr.to_csv(expr_out_path)

    # cnv
    ccle_cnv_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/ccle.lung_CNV_allgene.csv')
    utsw_cnv_path = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/utsw.lung_CNV_allgene.csv')
    cnv_out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_CNV_allgene.csv')

    ccle_cnv = pd.read_csv(ccle_cnv_path, index_col=0)
    utsw_cnv = pd.read_csv(utsw_cnv_path, index_col=0)

    cnv = itg.mergeExpression(CCLE=ccle_cnv, UTSW=utsw_cnv, keep_dup=['CCLE', 'UTSW'])
    # d = cnv.dropna().values.reshape(1,-1)
    # sns.distplot(d[(d >= -2) & (d <= 2)])
    # plt.savefig(os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/fig/ccle_utsw.cnv_dist_allgene.png'))

    cnv.to_csv(cnv_out_path)