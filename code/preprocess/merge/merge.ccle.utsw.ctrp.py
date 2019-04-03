##############################################################################################
###                               merge.ccle.utsw.ctrp.py                                  ###
##############################################################################################
# merge ccle and utsw lung cancer dataset

proj_path = '/work/bioinformatics/s418336/projects/DLMed'
import os
import re
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util
from utility import integrate as itg
from utility import plot as p
import matplotlib.pyplot as plt
import seaborn as sns

geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
# mutation
ccle_mut_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_Mutation_cancergene.csv')
utsw_mut_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_Mutation_cancergene.csv')
# expression
ccle_expr_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_RNAseq_cancergene.csv')
utsw_expr_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_RNAseq_cancergene.csv')
# cnv
ccle_cnv_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_CNV_cancergene.csv')
utsw_cnv_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_CNV_cancergene.csv')
# drug
ccle_drug_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung.gdsc.csv')
utsw_drug_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_Compound_ED50.csv')
ctrp_drug_path = os.path.join(proj_path, 'data/curated/Lung/ctrp.lung.csv')

dataset = 'ccle_utsw'

################################    main    ###################################
if dataset == 'ccle_utsw_ctrp':
    # >>>>>>>>>>>>>>>>   mutation   <<<<<<<<<<<<<<<<<<< #
    mut_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_Mutation_cancergene.csv')

    ccle_mutation = pd.read_csv(ccle_mut_path)
    utsw_mutation = pd.read_csv(utsw_mut_path)
    genes = pd.read_csv(geneanno_path)['Gene'].tolist()
    ccle_mutation = ccle_mutation.loc[ccle_mutation.Gene.isin(genes),:]
    utsw_mutation = utsw_mutation.loc[utsw_mutation.Gene.isin(genes),:]

    mutation = itg.mergeMutation(CCLE=ccle_mutation, UTSW=utsw_mutation)
    mutation.to_csv(mut_out_path, index=None)

    # >>>>>>>>>>>>>>>>   expression   <<<<<<<<<<<<<<<<<<< #
    expr_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_RNAseq_cancergene.csv')

    ccle_expr = pd.read_csv(ccle_expr_path, index_col=0)
    utsw_expr = pd.read_csv(utsw_expr_path, index_col=0)

    expr = itg.mergeExpression(CCLE=ccle_expr, UTSW=utsw_expr, keep_dup=['CCLE', 'UTSW'])
    expr.to_csv(expr_out_path)

    # >>>>>>>>>>>>>>>>   cnv   <<<<<<<<<<<<<<<<<<< #
    cnv_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_CNV_cancergene.csv')

    ccle_cnv = pd.read_csv(ccle_cnv_path, index_col=0)
    utsw_cnv = pd.read_csv(utsw_cnv_path, index_col=0)

    cnv = itg.mergeExpression(CCLE=ccle_cnv, UTSW=utsw_cnv, keep_dup=['CCLE', 'UTSW'])
    cnv.to_csv(cnv_out_path)

    # >>>>>>>>>>>>>>>>   drug data   <<<<<<<<<<<<<<<<<<< #
    drug_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_drug.csv')
    dens_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_density.png')

    ccle_drug = pd.read_csv(ccle_drug_path, usecols=['Aliases', 'Drug', 'IC50'])
    utsw_drug = pd.read_csv(utsw_drug_path)
    ctrp_drug = pd.read_csv(ctrp_drug_path)
    max_value = ccle_drug.IC50.max()
    min_value = ccle_drug.IC50.min()
    print(max_value, min_value)
    ccle_drug = util.cleanDrugData(ccle_drug)
    utsw_drug = util.cleanDrugData(utsw_drug, max_value=50, min_value=1e-4)
    ctrp_drug = util.cleanDrugData(ctrp_drug, max_value=max_value, min_value=1e-4)
    # merge
    drug = itg.mergeDrug(CCLE=ccle_drug, UTSW=utsw_drug, CTRP=ctrp_drug, keep_dup=['CCLE', 'UTSW', 'CTRP'])
    # plot
    f, ax = plt.subplots(figsize=(12,9))
    p.plot_density(CCLE=ccle_drug.LOG50, UTSW=utsw_drug.LOG50, CTRP=ctrp_drug.LOG50, MERGED=drug.LOG50, title='Density plot of drug log(IC50)', legend_title='Source', ax=ax)
    f.savefig(dens_path)
    drug.to_csv(drug_out_path, index=False)

elif dataset == 'ccle_utsw':
    # >>>>>>>>>>>>>>>>   mutation   <<<<<<<<<<<<<<<<<<< #
    mut_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_Mutation_cancergene.array.csv')

    ccle_mutation = pd.read_csv(ccle_mut_path)
    utsw_mutation = pd.read_csv(utsw_mut_path)
    genes = pd.read_csv(geneanno_path)['Gene'].tolist()
    ccle_mutation = ccle_mutation.loc[ccle_mutation.Gene.isin(genes),:]
    utsw_mutation = utsw_mutation.loc[utsw_mutation.Gene.isin(genes),:]

    mutation = itg.mergeMutation(CCLE=ccle_mutation, UTSW=utsw_mutation)[['Cell', 'Gene']].drop_duplicates()
    mutation['value'] = 1
    mutation = mutation.pivot(index='Cell', columns='Gene', values='value').fillna(0)
    mutation.index.name = mutation.columns.name = None
    mutation.to_csv(mut_out_path)

    # >>>>>>>>>>>>>>>>   expression   <<<<<<<<<<<<<<<<<<< #
    expr_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_RNAseq_cancergene.csv')

    ccle_expr = pd.read_csv(ccle_expr_path, index_col=0)
    utsw_expr = pd.read_csv(utsw_expr_path, index_col=0)

    expr = itg.mergeExpression(CCLE=ccle_expr, UTSW=utsw_expr, keep_dup=['CCLE', 'UTSW'])
    expr.to_csv(expr_out_path)

    # >>>>>>>>>>>>>>>>   cnv   <<<<<<<<<<<<<<<<<<< #
    cnv_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/merged.lung_CNV_cancergene.csv')

    ccle_cnv = pd.read_csv(ccle_cnv_path, index_col=0)
    utsw_cnv = pd.read_csv(utsw_cnv_path, index_col=0)

    cnv = itg.mergeExpression(CCLE=ccle_cnv, UTSW=utsw_cnv, keep_dup=['CCLE', 'UTSW'])
    cnv.to_csv(cnv_out_path)

    # >>>>>>>>>>>>>>>>   drug data   <<<<<<<<<<<<<<<<<<< #
    drug_out_path = os.path.join(proj_path, 'data/curated/Lung/merged/ccle_utsw.lung_drug.csv')
    dens_path = os.path.join(proj_path, 'data/curated/Lung/merged/ccle_utsw.lung_density.png')

    ccle_drug = pd.read_csv(ccle_drug_path, usecols=['Aliases', 'Drug', 'IC50'])
    utsw_drug = pd.read_csv(utsw_drug_path)
    max_value = ccle_drug.IC50.max()
    min_value = ccle_drug.IC50.min()
    print(max_value, min_value)
    ccle_drug = util.cleanDrugData(ccle_drug)
    utsw_drug = util.cleanDrugData(utsw_drug, max_value=50, min_value=min_value)
    # merge
    drug = itg.mergeDrug(CCLE=ccle_drug, UTSW=utsw_drug, keep_dup=['CCLE', 'UTSW'])
    # plot
    f, ax = plt.subplots(figsize=(12,9))
    p.plot_density(CCLE=ccle_drug.LOG50, UTSW=utsw_drug.LOG50, MERGED=drug.LOG50, title='Density plot of drug log(IC50)', legend_title='Source', ax=ax)
    f.savefig(dens_path)
    drug.to_csv(drug_out_path, index=False)
