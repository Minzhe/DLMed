################################################################################
###                           curate.compound.py                             ###
################################################################################

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

compound_path = os.path.join(proj_path, 'data/UTSW_MW/Cell.Line.Compound.xlsx')
out_path = os.path.join(proj_path, 'data/curated/Lung/utsw.mw/Lung_Compound_ED50.csv')

############################    main    ##############################
drug_ED50 = pd.read_excel(compound_path, sheet_name='ED50')
drug_ED50 = drug_ED50.iloc[:,np.concatenate([[0], list(range(7, drug_ED50.shape[1]))])]
drug_ED50 = pd.melt(drug_ED50, id_vars=['SW_ID'], var_name='Cell', value_name='ED50')
drug_ED50.Cell = drug_ED50.Cell.apply(util.cleanCellLine)
drug_ED50 = drug_ED50[['Cell', 'SW_ID', 'ED50']]
drug_ED50.columns = ['Cell', 'Compound', 'ED50']
drug_ED50.sort_values(by=['Cell', 'Compound'])
drug_ED50 = drug_ED50.loc[drug_ED50.ED50.notnull(),:]
drug_ED50.to_csv(out_path, index=None)