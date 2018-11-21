#########################################################################################
###                                  curate.expr.py                                   ###
#########################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os
import pandas as pd
import numpy as np

luad_path = os.path.join(proj_dir, 'data/TCGA/LUAD.rnaseqv2_unc_Level_3__RSEM_genes_normalized.txt')
lusc_path = os.path.join(proj_dir, 'data/TCGA/LUSC.rnaseqv2_unc_Level_3__RSEM_genes_normalized.txt')
expr_out = os.path.join(proj_dir, 'data/curated/Lung/tcga/Lung.expr.csv')
geneanno_path = os.path.join(proj_dir, 'data/COSMIC/Census_all.csv')

###############################  main  ###################################
gene_list = list(pd.read_csv(geneanno_path, usecols=[0], squeeze=True))

luad_expr = pd.read_table(luad_path, sep='\t', index_col=0, skiprows=[1])
luad_expr.index.name = None
luad_expr.index = luad_expr.index.to_series().apply(lambda x: x.split('|')[0])
luad_expr = luad_expr.loc[luad_expr.index.to_series().isin(gene_list),:]
luad_expr = luad_expr.T
luad_expr = np.log2(luad_expr + 1)
luad_expr['Type'] = 'LUAD'

lusc_expr = pd.read_table(lusc_path, sep='\t', index_col=0, skiprows=[1])
lusc_expr.index.name = None
lusc_expr.index = lusc_expr.index.to_series().apply(lambda x: x.split('|')[0])
lusc_expr = lusc_expr.loc[lusc_expr.index.to_series().isin(gene_list),:]
lusc_expr = lusc_expr.T
lusc_expr = np.log2(lusc_expr + 1)
lusc_expr['Type'] = 'LUSC'

expr = pd.concat([luad_expr, lusc_expr], axis=0)
expr.to_csv(expr_out)

