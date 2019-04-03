#################################################################################
###                              curate.CNV.py                                ###
#################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import pandas as pd

data_path = os.path.join(proj_dir, 'data/PDX/CNV.csv')
out_path = os.path.join(proj_dir, 'data/curated/PDX/PDX.CNV.allgenes.csv')
gene_set = 'all_gene'

###############################   main   #################################
expr = pd.read_csv(data_path, index_col=0).T
expr.to_csv(out_path)