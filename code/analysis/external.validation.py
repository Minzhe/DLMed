####################################################################################################################
###                                         external.validation.py                                               ###
####################################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
# proj_dir = 'Z:/bioinformatics/s418336/projects/DLMed'
import os
import numpy as np
import sys
import re
sys.path.append(os.path.join(proj_dir, 'code'))
import pandas as pd
import utility.utility as util
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

#################################    main   #####################################
def read_john_avail(path):
    john = pd.read_csv(path, usecols=[1,3]+list(range(41,124)))
    john = john.loc[john['Tumor Type'] == 'Lung',:].drop(['Tumor Type'], axis=1)
    john.columns = john.columns.to_series().apply(lambda x: re.sub(r'[-,\']', '', x.split(' ')[0]).upper())
    john['CELL'] = john['CELL'].apply(util.cleanCellLine)
    john = john.melt(id_vars=['CELL'], var_name='Drug', value_name='Source')
    john = john.loc[john.Source.notnull(),:]
    john.columns = ['Cell', 'Drug', 'Source']
    john.Source = 'LUKE'
    return john

def read_john_ic50(path):
    john = pd.read_csv(path, usecols=[1,3]+list(range(41,124)))
    john = john.loc[john['Tumor Type'] == 'Lung',:].drop(['Tumor Type'], axis=1)
    # convert unit
    for col in list(john.columns[1:]):
        if '(uM)' in col:
            del john[col]
        elif '(nM)' in col:
            john[col] = john[col] - 3
        elif 'mM' in col:
            del john[col]
        elif '(ug/ml)' in col:
            del john[col]
        else:
            raise ValueError('Unrecognized column {}'.format(col))
    john.columns = john.columns.to_series().apply(lambda x: re.sub(r'[-,\']', '', x.split(' ')[0]).upper())
    john['CELL'] = john['CELL'].apply(util.cleanCellLine)
    john = john.melt(id_vars=['CELL'], var_name='Drug', value_name='Luke_LOGIC50')
    john = john.loc[john.Luke_LOGIC50.notnull(),:]
    john['Luke_LOGIC50'] = john['Luke_LOGIC50'] * np.log(10)
    return john


def merge_gdsc_john(john, gdsc, keep='both'):
    gdsc = gdsc.loc[gdsc['Train_split_point'].notnull(),['Cell', 'Drug', 'Source']]
    gdsc['Source'][gdsc['Source'] == 'CCLE'] = 'GDSC'
    if keep == 'both':
        john = john.loc[john.Cell.isin(gdsc.Cell) & john.Drug.isin(gdsc.Drug)]
    elif keep != 'all':
        raise ValueError('Keep should be both or all.')
    data = pd.concat([gdsc, john])
    data = data.pivot_table(values='Source', index='Cell', columns='Drug', aggfunc=lambda x: ','.join(x))
    # print((gdsc == 'LUKE').values.sum()) # 886
    return data

def plot_correlation_val(pred, truth, out_dir):
    pred.columns = ['Cell', 'Drug', 'Truth', 'Pred']
    truth.columns = ['Cell', 'Drug', 'Truth']
    truth = truth.loc[truth.Cell.isin(pred.Cell) & truth.Drug.isin(pred.Drug),:]
    pred = pred.loc[pred.Truth.isnull(),:].drop(['Truth'], axis=1)
    # merge
    data = pred.merge(truth, on=['Cell', 'Drug'], how='inner')
    data_path = os.path.join(out_dir, 'validation_nm.csv')
    img_path = os.path.join(out_dir, 'validation_nm.png')
    data.to_csv(data_path, index=None)
    # plot
    lims = [data[['Pred','Truth']].values.min()-0.5, data[['Pred','Truth']].values.max()+0.5]
    plt.subplots(figsize=(12,12))
    plt.scatter(data['Pred'], data['Truth'], s=5)
    plt.title('Pred vs LUKE')
    plt.xlabel('Predicted log(ic50)')
    plt.ylabel('Observed log(ic50)')
    plt.plot(lims, lims, 'r--', alpha=0.75, zorder=0)
    plt.savefig(img_path)

def plot_correlation_true(pred, truth, out_dir):
    pred.columns = ['Cell', 'Drug', 'CCLE', 'Pred']
    truth.columns = ['Cell', 'Drug', 'LUKE']
    truth = truth.loc[truth.Cell.isin(pred.Cell) & truth.Drug.isin(pred.Drug),:]
    pred = pred.loc[pred.CCLE.notnull(),:].drop(['Pred'], axis=1)
    # merge
    data = pred.merge(truth, on=['Cell', 'Drug'], how='inner')
    data_path = os.path.join(out_dir, 'validation_true_nm.csv')
    data.to_csv(data_path)
    img_path = os.path.join(out_dir, 'validation_true_nm.png')
    # plot
    lims = [data[['CCLE','LUKE']].values.min()-0.5, data[['CCLE','LUKE']].values.max()+0.5]
    plt.subplots(figsize=(12,12))
    plt.scatter(data['CCLE'], data['LUKE'], s=5)
    plt.title('CCLE vs LUKE')
    plt.xlabel('CCLE')
    plt.ylabel('LUKE')
    plt.plot(lims, lims, 'r--', alpha=0.75, zorder=0)
    plt.savefig(img_path)

def plot_correlation_sub(df, drop_cell, drop_drug, title, out_dir):
    if drop_cell is not None:
        df = df.loc[~df['Cell'].isin(drop_cell),:]
    if drop_drug is not None:
        df = df.loc[~df['Drug'].isin(drop_drug),:]
    # plot
    r2 = round(r2_score(df['Truth'], df['Pred']), 3)
    lims = df[['Pred', 'Truth']].values.min()-0.5, df[['Pred', 'Truth']].values.max()+0.5
    lims = [-9, 6]
    plt.subplots(figsize=(12,12))
    plt.scatter(df['Pred'], df['Truth'], s=20)
    plt.title(title); plt.xlabel('Predicted log.ic50'); plt.ylabel('Observed log.ic50')
    plt.plot(lims, lims, 'r--', alpha=0.75, zorder=0)
    # plt.text(lims[0], lims[1], s='r2_score: {}'.format(str(r2)))
    plt.grid(linestyle='--')
    plt.savefig(os.path.join(out_dir, 'img/{}.png'.format(title)))
    df.to_csv(os.path.join(out_dir, '{}.csv'.format(title)), index=None)


def seprateCorr(df, col, val_dir, cutoff, write):
    mse = df.groupby(by=[col]).apply(lambda x: calMSE(x))
    corr = df.groupby(by=[col]).apply(lambda x: calR2(x))
    res = pd.concat([mse, corr], axis=1)
    res.columns = ['MSE', 'R2']
    if write:
        res.to_csv(os.path.join(val_dir, 'corr.{}.csv'.format(col)))
    bad = res.index[res['R2'] < cutoff]
    return bad

def calMSE(df):
    return np.mean(np.square(df['Truth'] - df['Pred']))

def calR2(df):
    return r2_score(df['Truth'], df['Pred'])

#################################    main   #####################################
'''
### clinical information
infer_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/cell_inference.csv')
utsw_path = os.path.join(proj_dir, 'data/UTSW_MW/Cell.Line.Database v1.73.csv')
out_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/inference.clinical.csv')
infer_cell = pd.read_csv(infer_path, squeeze=True)
cols = ['Cell Line', 'Tumor Type', 'Tumor Subtype', 'Reclassification using Expression Signature', 'EMT Phenotype', 'Age', 'Race', 'Gender', 'Tumor Source', 
        'Anatomical Site', 'Notes', 'SMOKER', 'SMOKINGPYR', 'STAGE', 'Date Cell Line', 'DATE1STRX', 'PRIOR_RX', 'CT', 'RT', 'RX_AFTER', 'DRUGS', 'RESP_BEST',
        'PROTOCOL', 'LRXDATE', 'LDRUGS', 'BEST_RESP', 'DEATH', 'DXDATE', 'CELLDDC', 'SCSUBTYPE', ]
utsw = pd.read_csv(utsw_path, usecols=cols)
utsw['Cell Line'] = utsw['Cell Line'].apply(util.cleanCellLine)
utsw = utsw[utsw['Cell Line'].isin(infer_cell)].sort_values(by=['Cell Line'])
utsw.to_csv(out_path, index=False)
'''

### john lab data
model = 'full_cnn'
john_path = os.path.join(proj_dir, 'data/UTSW_MW/Cell.Line.Database v1.73.csv')
gdsc_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv')
pred_path = os.path.join(proj_dir, 'result/prediction/{}/results.csv'.format(model))
out_path  = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/data.availability_validation.csv')
val_dir   = os.path.join(proj_dir, 'result/prediction/{}'.format(model))

# merge john gdsc data
# john = read_john_avail(john_path)
# gdsc = pd.read_csv(gdsc_path)
# merge_gdsc_john(john, gdsc)

# plot correlation
# john = read_john_ic50(john_path)
# pred = pd.read_csv(pred_path)
# plot_correlation_val(pred, john, val_path)
# plot_correlation_true(pred, john, val_path)

### plot sub correlation
model = 'full_cnn'
val_path = os.path.join(proj_dir, 'result/prediction/{}/validation.csv'.format(model))
val = pd.read_csv(val_path)
val = val.loc[val['Truth'] != 0,:]
bad_cell = seprateCorr(val, 'Cell', val_dir, cutoff=0.7, write=False)
# bad_drug = seprateCorr(val, 'Drug', val_dir, cutoff=-6, write=False)
plot_correlation_sub(val, bad_cell, None, 'Pred vs Truth', val_dir)