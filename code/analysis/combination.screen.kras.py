#########################################################################################################
###                                   combination.effect.kras.py                                      ###
#########################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import re
import glob
import pandas as pd
import numpy as np
import statsmodels.api as sm
import multiprocessing as mp

####################################   function  ####################################
# *************************************************** #
#               parse prediction result               #
# *************************************************** #
def parseAllCombination(dir, out_dir, n_job):
    paths = glob.glob(dir)
    n = len(paths)
    print('Number of files: {}'.format(n))
    for i in range(0, n, n_job):
        tmp_paths = paths[i:i+n_job]
        jobs = []
        for path in tmp_paths:
            cell = re.findall(r'cell_([0-9]+)_drug', path)[0]
            out8_path = os.path.join(out_dir, 'combination.kras_8.cell_{}.csv'.format(cell))
            out9_path = os.path.join(out_dir, 'combination.kras_9.cell_{}.csv'.format(cell))
            p = mp.Process(target=parseCombination, args=(path, out8_path, out9_path,))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()

def parseCombination(path, out8_path, out9_path):
    print('Reading data {}'.format(path.split('/')[-1]))
    data = pd.read_csv(path)
    data['log_ic50'] = round(data['log_ic50'], 6)
    # analyze 359_8 mutation
    print('Analyzing Kras_8 mutation ...')
    data8 = data.loc[data['359_9'] == 0,:].drop(['359_9'], axis=1)
    data8 = reduceCombTable(data8)
    data8.to_csv(out8_path, index=0)
    # analyze 359_9 mutation
    print('Analyzing Kras_9 mutation ...')
    data9 = data.loc[data['359_8'] == 0,:].drop(['359_8'], axis=1)
    data9 = reduceCombTable(data9)
    data9.to_csv(out9_path, index=0)

def reduceCombTable(df):
    df.columns = ['kras'] + list(df.columns)[1:]
    df['cond'] = df['kras'] + df['mut']*2
    df = df.drop(['kras', 'mut'], axis=1)
    df = df.pivot_table(values='log_ic50', index=['cellline_id', 'drug_id', 'row', 'col'], columns='cond').reset_index()
    df.columns = ['cell_id', 'drug_id', 'gene_id', 'win', 'no_mut', 'mut_kras', 'mut_gene', 'mut_both']
    df.columns.name = None
    return df

# ***************************************************** #
#               aggregate location to gene              #
# ***************************************************** #
def aggregateAllCombination(path, out_path, n_jobs):
    paths = glob.glob(path)
    data = list()
    for i in range(0, len(paths), n_jobs):
        tmp_data = mp.Manager().dict()
        jobs = []
        print('.........')
        for j in range(i, min(i+n_jobs, len(paths))):
            tmp_path = paths[j]
            print('Analyzing {}'.format(tmp_path.split('/')[-1]))
            p = mp.Process(target=aggregateCombination, args=(tmp_path, tmp_data,))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()
        tmp_data = pd.concat(list(tmp_data.values()))
        data.append(tmp_data)
    # merge
    data = pd.concat(data)
    data = data.groupby(by=['cell_id', 'drug_id', 'gene_id'], as_index=False).agg(lambda x: max(x, key=abs))
    data.to_csv(out_path, index=False)

def aggregateCombination(path, return_dict):
    data = pd.read_csv(path)
    data['mut_kras'] = data['mut_kras'] - data['no_mut']
    data['mut_gene'] = data['mut_gene'] - data['no_mut']
    data['mut_both'] = data['mut_both'] - data['no_mut']
    data['no_mut'] = data['no_mut'] - data['no_mut']
    data = data.drop(['win'], axis=1)
    data = data.groupby(by=['cell_id', 'drug_id', 'gene_id'], as_index=False).agg(lambda x: max(x, key=abs))
    return_dict[path] = data

# ******************************************************** #
#               calculate combination effect               #
# ******************************************************** #
def calCombinationEffect(path, out_path, l_proc, l_chunk):
    paths = glob.glob(path)
    n = sum(1 for row in open(paths[1], 'r'))
    n_split = -(-n // l_chunk)
    for i in range(n_split):
        # multiprocessing
        n_jobs = -(-l_chunk // l_proc)
        manager = mp.Manager()
        return_dict = manager.dict()
        jobs = []
        for j in range(n_jobs):
            skip = max(1, j*l_proc+i*l_chunk)
            p = mp.Process(target=calCombinationEffect_process, args=(paths, skip, l_proc, return_dict))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()
        comb_eff = pd.concat(list(return_dict.values()), axis=0)
        tmp_out = out_path.strip('.csv') + '.{}.csv'.format(i)
        comb_eff.to_csv(tmp_out, index=None)
                
def calCombinationEffect_process(paths, skip, lines, return_dict):
    print('Analyzing line {} to {}'.format(skip, skip+lines))
    data = list()
    for path in paths:
        tmp = pd.read_csv(path, header=None, skiprows=skip, nrows=lines, names=['cell_id', 'drug_id', 'gene_id', 'win', 'no_mut', 'mut_kras', 'mut_gene', 'mut_both'])
        data.append(tmp)
    data = pd.concat(data)
    data.index = range(data.shape[0])
    data['mut_kras'] = data['mut_kras'] - data['no_mut']
    data['mut_gene'] = data['mut_gene'] - data['no_mut']
    data['mut_both'] = data['mut_both'] - data['no_mut']
    data['no_mut'] = data['no_mut'] - data['no_mut']
    comb_eff = list()
    # calculate combination effect
    comb_eff = data.groupby(by=['drug_id', 'gene_id', 'win']).apply(lambda x: fitLinearModel(x[['no_mut', 'mut_kras', 'mut_gene', 'mut_both']]))
    comb_eff.columns = ['coef', 'neg.log.p']
    comb_eff = comb_eff.reset_index()
    return_dict[(skip, lines)] = comb_eff

def calCombinationEffect_func(df):
    print('Calculating combination effect ...')
    comb_eff = df.groupby(by=['drug_id', 'gene_id']).apply(lambda x: fitLinearModel(x[['no_mut', 'mut_kras', 'mut_gene', 'mut_both']]))
    comb_eff.columns = ['coef', 'neg.log.p']
    comb_eff = comb_eff.reset_index()
    return comb_eff

def fitLinearModel(df):
    n = df.shape[0]
    x1 = np.array([0]*n + [1]*n + [0]*n + [1]*n)
    x2 = np.array([0]*n + [0]*n + [1]*n + [1]*n)
    x12 = np.array([0]*n + [0]*n + [0]*n + [1]*n)
    X = np.array([x1, x2, x12]).T
    y = np.array(df).T.reshape(-1)
    ols = sm.OLS(y, X).fit()
    return pd.Series([round(ols.params[2], 6), round(-np.log(ols.pvalues[2]), 6)])

def getDrugGeneWin(path):
    data = pd.read_csv(path, usecols=['drug_id', 'gene_id', 'win'])
    return (zip(data['drug_id'], data['gene_id'], data['win']))

####################################    main   ########################################
model = 'full_cnn'
if model == 'gene_cnn':
    simu_path = '/project/bioinformatics/Xiao_lab/shared/Zhan_Xu/CCLE_simulation/res_simu_supp_on_cellline_mask_kras/*.csv'
    out_dir = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.gene_cnn')

    # parse combination
    parseAllCombination(simu_path, out_dir, n_job=16)
    # aggregation loci to gene
    kras_path = os.path.join(out_dir, 'combination.kras_*.cell_*.csv')
    out_path = os.path.join(out_dir, 'combination.kras.agg.csv')
    aggregateAllCombination(kras_path, out_path, n_jobs=8)
    # calculate combination effect
    eff_path = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.gene_cnn.analysis/kras.combination.eff.csv')
    data = pd.read_csv(out_path)
    data = calCombinationEffect_func(data)
    data.to_csv(eff_path, index=False)

elif model == 'full_cnn':
    simu_path = '/project/bioinformatics/Xiao_lab/shared/Zhan_Xu/CCLE_simulation/res_simu_on_cellline_mask_kras/*.csv'
    out_dir = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.full_cnn')

    # # parse combination
    # parseAllCombination(simu_path, out_dir, n_job=16)
    # # aggregation loci to gene
    # kras_path = os.path.join(out_dir, 'combination.kras_*.cell_*.csv')
    out_path = os.path.join(out_dir, 'combination.kras.agg.csv')
    # aggregateAllCombination(kras_path, out_path, n_jobs=16)
    # calculate combination effect
    eff_path = os.path.join(proj_dir, 'result/simulation/simulation.combination_screen_kras.full_cnn.analysis/kras.combination.eff.csv')
    data = pd.read_csv(out_path)
    data = calCombinationEffect_func(data)
    data.to_csv(eff_path, index=False)


# kras8_path = os.path.join(out_dir, 'combination.kras_8.*.csv')
# kras9_path = os.path.join(out_dir, 'combination.kras_9.*.csv')
# out8_path = os.path.join(proj_dir, 'result/simulation/kras_combination_simulation_analysis/kras8_combination_effect.csv')
# calCombinationEffect(kras8_path, out8_path, l_proc=10000, l_chunk=10000*16)
# out9_path = os.path.join(proj_dir, 'result/simulation/kras_combination_simulation_analysis/kras9_combination_effect.csv')
# calCombinationEffect(kras9_path, out9_path, l_proc=10000, l_chunk=10000*10)
