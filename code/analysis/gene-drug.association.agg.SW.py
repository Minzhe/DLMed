####################################################################################################
###                                   gene-drug.association.py                                   ###
####################################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import glob
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF

###############################   function  ##################################
def agg_gene(df, agg):
    data = df.copy()
    data.index = data.index.to_series().apply(lambda x: x.split('_')[0])
    if agg == 'max':
        data = data.groupby(by=data.index).agg(lambda x: max(x, key=abs))
    elif agg == 'median':
        data = data.groupby(by=data.index).agg(np.median)
    else:
        raise ValueError('Unrecognized agg function.')
    return data

def calculate_p(df, drug):
    sw_df = df[drug]
    df = abs(df)
    # empirical distribution
    p_val = list()
    for idx, row in df.iterrows():
        ecdf = ECDF(list(row))
        p = 1 - ecdf(abs(sw_df.loc[idx,:]))
        p_val.append(p)
    p_val = pd.DataFrame(p_val, index=sw_df.index, columns=[item+'_p' for item in drug])
    sw_df = round(pd.concat([sw_df, p_val], axis=1).sort_index(axis=1), 6)
    return sw_df

def summarize_drug(drugs, out_dir, **data):
    for drug in drugs:
        out_path = os.path.join(out_dir, 'gene-drug.association.{}.csv'.format(drug))
        drug_df = []
        for name in data.keys():
            df = data[name]
            colnames = df.columns.tolist()
            colnames = [col for col in colnames if drug in col]
            df = df[colnames]
            df.columns = [col+'_'+name for col in colnames]
            drug_df.append(df)
        drug_df = pd.concat(drug_df, axis=1)
        drug_df = drug_df.loc[drug_df.iloc[:,1] < 0.05,:]
        drug_df = drug_df.loc[drug_df.iloc[:,1].sort_values().index,:]
        drug_df.to_csv(out_path)
        



###############################   main  ##################################
diff_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/random_single_mut_simulation.gene-drug.association_median.csv')
# out_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis/random_single_mut_simulation.gene-drug.association_median.SW_drug.csv')
out_dir = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis')
drug = ['SW036310', 'SW022906', 'SW069087']

diff = pd.read_csv(diff_path, index_col=0)
diff_max = agg_gene(diff, agg='max')
diff_median = agg_gene(diff, agg='median')
diff_max = calculate_p(diff_max, drug)
diff_median = calculate_p(diff_median, drug)
summarize_drug(drug, out_dir, max=diff_max, median=diff_median)
# sw_diff = sw_diff.loc[((diff['SW036310'] > 0.1) | (diff['SW022906'] > 0.1) | (diff['SW069087'] > 0.1)) & ((sw_diff['SW036310_p'] < 0.05) | (sw_diff['SW022906_p'] < 0.05) | (sw_diff['SW069087_p'] < 0.05)),:]
# sw_diff.to_csv(out_path)

##############   raw cell line data  ################
# file_path = '/project/bioinformatics/Xiao_lab/shared/DL_drug_prediction/single_mutation_diff_logic50/*.csv'
# file_path = glob.glob(file_path)
# sw036310, sw022906, sw069087= dict(), dict(), dict()
# for path in file_path:
#     print(path)
#     cell = path.split('.')[0].split('_')[-1]
#     data = pd.read_csv(path, index_col=0)
#     sw036310[cell] = list(data['SW036310'])
#     sw022906[cell] = list(data['SW022906'])
#     sw069087[cell] = list(data['SW069087'])
# sw036310 = pd.DataFrame(sw036310, index=data.index)
# sw022906 = pd.DataFrame(sw022906, index=data.index)
# sw069087 = pd.DataFrame(sw069087, index=data.index)
# sw036310.to_csv(os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW036310.csv'))
# sw022906.to_csv(os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW022906.csv'))
# sw069087.to_csv(os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis_SW/SW069087.csv'))