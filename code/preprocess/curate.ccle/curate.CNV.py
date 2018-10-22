###########################################################################
###                           curate.CNV.py                             ###
###########################################################################

proj_path = 'D:/projects/DLCell'
import numpy as np
import pandas as pd
import os
import re
import sys
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util
from utility import plot as p
import matplotlib.pyplot as plt

cnv_path = os.path.join(proj_path, "data/DepMap/DepMap.CNV.2018.q3.csv")
geneanno_path = os.path.join(proj_path, "data/curated/cancer.gene.anno.csv")
cellline_path = os.path.join(proj_path, "data/curated/lung/ccle/cellline.meta.ccle.lung.csv")
out_path = os.path.join(proj_path, "data/curated/lung/ccle/lung_CNV_cancergene.csv")

############################   main   ################################
# read cancer gene list
anno = pd.read_csv(geneanno_path)['Gene'].tolist()

# read cell line info
cellline = pd.read_csv(cellline_path)
geneNameEncoder = util.geneEncoder(broadid=cellline['Broad_ID'].tolist(), alias=cellline['Aliases'].tolist())

# read cnv data
cnv = pd.read_csv(cnv_path, index_col=0)
cnv.columns = [gene.split('(')[0].strip() for gene in cnv.columns.values]
cnv = cnv[sorted(cnv.columns.tolist())]

# filter gene
# cnv = cnv.loc[:,cnv.columns.isin(anno)]
# print(anno)
cnv = cnv.loc[:,cnv.columns.to_series().isin(anno)]

# filter cell
cnv.index = geneNameEncoder.broad2name(cnv.index)
cnv = cnv.loc[cnv.index != 'NA',:]
cnv.sort_index(inplace=True)

# plot density
p.plot_density(cnv, xlim=(-4,4))
plt.show()

cnv.to_csv(out_path)

