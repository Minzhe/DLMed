###########################################################################
###                           curate.CNV.py                             ###
###########################################################################

proj_path = '/work/bioinformatics/s418336/projects/DLMed/'
import numpy as np
import pandas as pd
import os
import re
import sys
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util
from utility import plot as p
import matplotlib.pyplot as plt

genelist = 'allgene'

############################    main    ############################
if genelist == 'cancergene':
    cnv_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.CNV.csv')
    geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
    out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_CNV_cancergene.csv')
    # read cancer gene list
    cancer_genes = pd.read_csv(geneanno_path, usecols=['Gene'], squeeze=True)

    # read cnv data
    cnv = pd.read_csv(cnv_path)
    del cnv['Chr'], cnv['MapInfo'], cnv['Unigene'], cnv['Cytoband']
    cnv = cnv.loc[cnv.Symbol.isin(cancer_genes),:]
    cnv.set_index('Symbol', inplace=True)
    cnv.index.name = None
    cnv.sort_index(inplace=True)
    cnv = cnv.transpose()
    cnv.index = cnv.index.to_series().apply(util.cleanCellLine)
    cnv = cnv.groupby(by=[cnv.index]).mean()
    cnv.sort_index(inplace=True)
    cnv = cnv - 2

    cnv.to_csv(out_path)

############################    main    ############################
elif genelist == 'allgene':
    cnv_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.CNV.csv')
    out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/utsw.lung_CNV_allgene.csv')

    # read cnv data
    cnv = pd.read_csv(cnv_path)
    del cnv['Chr'], cnv['MapInfo'], cnv['Unigene'], cnv['Cytoband']
    cnv.set_index('Symbol', inplace=True)
    cnv.index.name = None
    cnv.sort_index(inplace=True)
    cnv = cnv.transpose()
    cnv.index = cnv.index.to_series().apply(util.cleanCellLine)
    cnv = cnv.groupby(by=[cnv.index]).mean()
    cnv.sort_index(inplace=True)
    cnv = cnv - 2

    cnv.to_csv(out_path)
