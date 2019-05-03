###########################################################################
###                           curate.DRC.py                             ###
###########################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import pandas as pd
import numpy as np
import utility.utility as util

###############################   function  ###############################
def readDRC(path, out_dir):
    with open(path, 'r') as f:
        lines = f.readlines()
        groups = groupLines(lines)
    for pair in groups:
        info = parsePair(pair)
        # print(info['drc'])
        # exit()
        out_path = os.path.join(out_dir, 'drc.{}-{}-{}.csv'.format(info['Cell'], info['Drug'], info['Exp']))
        write_output(info, out_path)
        
def parsePair(lines):
    pair = dict()
    pair['Exp'] = lines[0].split('\t')[-1]
    pair['Cell'] = util.cleanCellLine(lines[1].split('\t')[-1])
    pair['Unit'] = lines[4].split('\t')[-1]
    pair['Drug'] = util.cleanDrugName(lines[3].split('\t')[-1])
    pair['IC50'] = float(lines[6].split('\t')[-1])
    pair['ED50'] = float(lines[8].split('\t')[-1])
    pair['MaxRes'] = float(lines[10].split('\t')[-1])
    pair['drc'] = [','.join(line.strip().split('\t')[1:]) for line in lines[17:26]]
    return pair

### >>>>>>>>>>>>>>>>>>>>  utility  <<<<<<<<<<<<<<<<<<<<< ###
def groupLines(lines):
    groups, pair = [], []
    for line in lines:
        line = line.strip()
        if len(line) > 0:
            pair.append(line)
        else:
            groups.append(pair)
            pair = []
    return groups

def write_output(info, path):
    print(path)
    with open(path, 'w') as f:
        print('# Exp ID: {}'.format(info['Exp']), file=f)
        print('# Cell: {}'.format(info['Cell']), file=f)
        print('# Drug: {}'.format(info['Drug']), file=f)
        print('# IC50: {}'.format(info['IC50']), file=f)
        print('# ED50: {}'.format(info['ED50']), file=f)
        print('# Unit: {}'.format(info['Unit']), file=f)
        print('# MaxRes: {}'.format(info['MaxRes']), file=f)
        for line in info['drc']:
            print(line, file=f)

def parse_table(lines):
    dose = [[np.nan, np.nan] + list(map(float, lines[0].split('\t')[2:-1]))]
    arr = np.array([list(map(float, l.split('\t')[1:])) for l in lines[1:]])
    arr = np.concatenate([arr[:,[0,-1]], arr[:,1:-1]], axis=1)
    arr = np.concatenate([np.array(dose), arr])
    return arr

###############################    main   #################################
data_path = os.path.join(proj_dir, 'data/UTSW_MW/DRC/Drug response curves Minna lab 4-23-19.txt')
out_dir = os.path.join(proj_dir, 'data/curated/Lung/utsw.mw/DRC')
readDRC(data_path, out_dir)