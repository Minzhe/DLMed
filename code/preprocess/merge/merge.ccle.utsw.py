#########################################################################################
###                               merge.ccle.utsw.py                                  ###
#########################################################################################
# merge ccle and utsw lung cancer dataset

proj_path = 'D:/projects/DLMed'
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
ccle_mut_path = os.path.join(proj_path, 'data/curated/lung/ccle/ccle.lung_Mutation_cancergene.csv')
utsw_mut_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/utsw.lung_Mutation_cancergene.csv')
mut_out_path = os.path.join(proj_path, 'data/curated/lung/merged/merged.lung_Mutation_cancergene.csv')
# expression
ccle_expr_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_RNAseq_cancergene.csv')
utsw_expr_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_RNAseq_cancergene.csv')
expr_out_path = os.path.join(proj_path, 'data/curated/lung/merged/merged.lung_RNAseq_cancergene.csv')
# cnv
ccle_cnv_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung_CNV_cancergene.csv')
utsw_cnv_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_CNV_cancergene.csv')
cnv_out_path = os.path.join(proj_path, 'data/curated/lung/merged/marged.lung_CNV_cancergene.csv')
# drug
ccle_drug_path = os.path.join(proj_path, 'data/curated/Lung/ccle/ccle.lung.gdsc.csv')
utsw_drug_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.Lung_Compound_ED50.csv')
ctrp_drug_path = os.path.join(proj_path, 'data/curated/Lung/ctrp.lung.csv')


################################    main    ###################################
# >>>>>>>>>>>>>>>>   mutation   <<<<<<<<<<<<<<<<<<< #
# ccle_mutation = pd.read_csv(ccle_mut_path)
# utsw_mutation = pd.read_csv(utsw_mut_path)
# genes = pd.read_csv(geneanno_path)['Gene'].tolist()
# ccle_mutation = ccle_mutation.loc[ccle_mutation.Gene.isin(genes),:]
# utsw_mutation = utsw_mutation.loc[utsw_mutation.Gene.isin(genes),:]

# mutation = itg.mergeMutation(CCLE=ccle_mutation, UTSW=utsw_mutation)
# mutation.to_csv(mut_out_path, index=None)

# # >>>>>>>>>>>>>>>>   expression   <<<<<<<<<<<<<<<<<<< #
# ccle_expr = pd.read_csv(ccle_expr_path, index_col=0)
# utsw_expr = pd.read_csv(utsw_expr_path, index_col=0)

# expr = itg.mergeExpression(CCLE=ccle_expr, UTSW=utsw_expr, merge_seq=['CCLE', 'UTSW'])
# expr.to_csv(expr_out_path)

# >>>>>>>>>>>>>>>>   cnv   <<<<<<<<<<<<<<<<<<< #
# ccle_cnv = pd.read_csv(ccle_cnv_path, index_col=0)
# utsw_cnv = pd.read_csv(utsw_cnv_path, index_col=0)

# cnv = itg.mergeExpression(CCLE=ccle_cnv, UTSW=utsw_cnv, merge_seq=['CCLE', 'UTSW'])
# cnv.to_csv(cnv_out_path)


# >>>>>>>>>>>>>>>>   drug data   <<<<<<<<<<<<<<<<<<< #
ccle_drug = pd.read_csv(ccle_drug_path, usecols=['Aliases', 'Drug', 'IC50'])
utsw_drug = pd.read_csv(utsw_drug_path)
ctrp_drug = pd.read_csv(ctrp_drug_path)
ccle_drug = util.cleanDrugData(ccle_drug)
utsw_drug = util.cleanDrugData(utsw_drug)
ctrp_drug = util.cleanDrugData(ctrp_drug)
# print(ccle_drug.head())
# print(utsw_drug.head())
# print(ctrp_drug.head())
# print(len(set(ccle_drug.Cell)))
# print(len(set(utsw_drug.Cell)))
# print(len(set(ctrp_drug.Cell)))
drug, dup_drug = itg.mergeDrug(CCLE=ccle_drug, UTSW=utsw_drug, CTRP=ctrp_drug)
print(drug)
print(dup_drug)
# p.plot_density(CCLE=ccle_drug.LOG50, UTSW=utsw_drug.LOG50, CTRP=ctrp_drug.LOG50)
# plt.show()
# print(set(utsw_drug.Cell) - set(ccle_drug.Aliases))
# print(set(ctrp_drug.Cell) - set(ccle_drug.Aliases))
# print(len(set(ccle_drug.Drug.str.upper())))
# print(len(set(ctrp_drug.Drug.str.upper())))
# print(len(set(ccle_drug.Drug.str.upper()) & set(ctrp_drug.Drug.str.upper())))
# print(set(ccle_drug.Drug.str.upper()) & set(ctrp_drug.Drug.str.upper()))