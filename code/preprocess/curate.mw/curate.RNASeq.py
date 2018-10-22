###########################################################################
###                         curate.RNASeq.py                            ###
###########################################################################

proj_path = 'D:/projects/DLCell'
import numpy as np
import pandas as pd
import os
import re
import sys
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util


rnaseq_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.RNASeq.txt')
geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/lung_RNAseq_cancergene.csv')


############################   main   ################################
# read cancer gene list
cancer_genes = pd.read_csv(geneanno_path, usecols=['Gene'], squeeze=True)

# read rnaseq data
rnaseq = pd.read_table(rnaseq_path, sep='\t')
rnaseq.drop(['ENTREZGENEID'], axis=1, inplace=True)
rnaseq = rnaseq.loc[rnaseq['SYMBOL'].isin(cancer_genes),:]
rnaseq.set_index('SYMBOL', inplace=True)
rnaseq.index.name = None
rnaseq.sort_index(inplace=True)
rnaseq = rnaseq.transpose()
rnaseq.index = rnaseq.index.to_series().apply(util.cleanCellLine)
rnaseq.sort_index(inplace=True)
rnaseq.to_csv(out_path)