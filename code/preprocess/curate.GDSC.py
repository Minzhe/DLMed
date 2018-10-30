###################################################################################
###                              curate.GDSC.py                                 ###
###################################################################################
# curate GDSC drug sensitivity dataset

proj_dir = "D:/projects/DLMed"
import os
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.join(proj_dir, 'code'))
from utility import utility as util

gdsc_path = os.path.join(proj_dir, 'data/Cosmic/GDSC.fitted_dose_response.v17.3.csv')
lungcell_path = os.path.join(proj_dir, 'data/curated/Lung/ccle/cellline.meta.ccle.lung.csv')
out_path = os.path.join(proj_dir, 'data/curated/Lung/ccle.lung.gdsc.csv')


#############################     main    ################################
gdsc_data = pd.read_csv(gdsc_path, usecols=['CELL_LINE_NAME', 'DRUG_NAME', 'PUTATIVE_TARGET', 'LN_IC50'])
gdsc_data.columns = ['Aliases', 'Drug', 'Target', 'IC50']
gdsc_data['Aliases'] = gdsc_data['Aliases'].apply(util.cleanAlias)
gdsc_data.sort_values(by=['Aliases', 'Drug'], inplace=True)

### lung ccle
lungcell_drugic50 = pd.read_csv(lungcell_path, usecols=['Broad_ID', 'Aliases', 'Primary Disease'])
lungcell_drugic50 = lungcell_drugic50.join(gdsc_data.set_index('Aliases'), on='Aliases', how='inner')
lungcell_drugic50.IC50 = np.exp(lungcell_drugic50.IC50)
lungcell_drugic50.to_csv(out_path, index=False)

	