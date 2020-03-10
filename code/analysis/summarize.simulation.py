########################################################################################
###                            summarize.simulation.py                               ###
########################################################################################

import os
import pandas as pd
import pickle as pkl
from glob import glob

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
data_path = os.path.join(proj_dir, 'result/simulation/source')
out_path = os.path.join(proj_dir, 'result/simulation/top_genes.csv')

##########################    function    ###############################
def read_simulation(path):
    f = open(path, 'rb')
    data = pkl.load(f)
    df = data['summary'][['gene_name', 'ic50']].iloc[:8,:]
    df['cell_min'], df['cell_max'], df['drug'] = data['cellline_min'], data['cellline_max'], data['drug']
    df = df[['cell_min', 'cell_max', 'drug', 'gene_name', 'ic50']]
    df.columns = ['cell_min', 'cell_max', 'drug', 'gene', 'logic50']
    return df

##########################  main  ############################
files = glob(os.path.join(data_path, '*.pkl'))
summary = [read_simulation(path) for path in files]
summary = pd.concat(summary)
summary.sort_values(by=['cell_min', 'cell_max', 'drug'], inplace=True)
summary.to_csv(out_path, index=None)