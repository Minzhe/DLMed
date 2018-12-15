###############################################################################
###                         curate.microarray.py                            ###
###############################################################################

proj_path = 'D:/projects/DLCell'
import numpy as np
import pandas as pd
import os
import re
import sys
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util

microarray_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.microarray.csv')
geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/lung_microarray_cancergene.csv')

############################   main   ################################
# read cancer gene list
cancer_genes = pd.read_csv(geneanno_path, usecols=['Gene'], squeeze=True)

# read rnaseq data
microarray = pd.read_csv(microarray_path).iloc[:,32:-1]
microarray = microarray.loc[microarray['Symbol'].isin(cancer_genes),:]
microarray = microarray.groupby(by=['Symbol'], sort=True).median()
microarray.index.name = None
microarray = microarray.transpose()
microarray.index = microarray.index.to_series().apply(util.cleanCellLine)
microarray.sort_index(inplace=True)
microarray.to_csv(out_path)