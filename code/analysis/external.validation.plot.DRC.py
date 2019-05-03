####################################################################################################################
###                                      external.validation.plot.DRC.py                                         ###
####################################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import glob
import pandas as pd
import numpy as np
from utility.utility import DRC

####################################     function    ###################################
def plotDRC(cell, drug, drc_dir, pred, out_dir):
    # read prediction
    pred = pred.loc[(pred['Cell'] == cell) & (pred['Drug'] == drug),:]
    pred_logic50, _ = float(pred['Pred']), float(pred['Truth'])
    # read drc
    paths = glob.glob(os.path.join(drc_dir, 'drc.{}-{}*.csv'.format(cell, drug)))
    for path in paths:
        ic50, unit, maxres = read_drc_param(path)
        dose, resp = read_drc_table(path, unit)
        # estimate ic50
        drc = DRC(dose=dose, resp=resp)
        if drc.fit:
            log_ic50 = np.log(drc.get_ic50())
            if abs(log_ic50 - pred_logic50) < 1.5:
                print('{} - {} | predicted: {}, estimated: {}'.format(cell, drug, pred_logic50, log_ic50))
                out_path = os.path.join(out_dir, os.path.basename(path).strip('.csv') + '.png')
                title = '{} - {}'.format(cell, drug)
                drc.plot_drc(pred_logic50, out_path, title)

### >>>>>>>>>>>>>>>>>  utility  <<<<<<<<<<<<<<<<<<< ###
def read_drc_table(path, unit):
    # read drc
    dose = pd.read_csv(path, skiprows=7, nrows=1, header=None).T.squeeze().values
    resp = pd.read_csv(path, skiprows=8, header=None).values
    # subtract background
    blank = resp[:,dose == 'Blank'].mean(axis=1).reshape((-1,1))
    resp = (resp - blank)[:,dose != 'Blank']
    dose = dose[dose != 'Blank']
    # scale to percentage
    base = resp[:,dose == 0].mean(axis=1).reshape((-1,1))
    resp = (resp / base)[:,dose != 0]
    dose = dose[dose != 0]
    # reshape array
    if unit == 'nanomolar':
        dose = dose / 1000
    elif unit != 'micromolar':
        raise ValueError('Unrecognizable unit.')
    resp = resp[~np.isnan(resp).any(axis=1),:]
    dose = np.repeat(dose, resp.shape[0])
    resp = resp.T.reshape(-1)
    return dose, resp

def read_drc_param(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    ic50 = float(lines[3].split(':')[-1].strip())
    ed50 = float(lines[4].split(':')[-1].strip())
    unit = lines[5].split(':')[-1].strip()
    maxres = float(lines[6].split(':')[-1].strip())
    return ic50, unit, maxres


#######################################     main    #######################################
drc_dir = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/DRC/')
pred_path = os.path.join(proj_dir, 'result/prediction/full_cnn/validation.csv')
out_dir = os.path.join(proj_dir, 'result/prediction/full_cnn/img.drc')

pred = pd.read_csv(pred_path)
pair = pred.loc[abs(pred['Pred'] - pred['Truth']) < 0.5,['Cell', 'Drug']]
for idx, row in pair.iterrows():
    plotDRC(row['Cell'], row['Drug'], drc_dir, pred, out_dir)
