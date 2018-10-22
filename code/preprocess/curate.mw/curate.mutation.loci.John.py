##################################################################################
###                         curate.mutation.loci.py                            ###
##################################################################################
# curate utsw mike white mutation in loci level

proj_path = 'D:/projects/DLCell'
import os
import re
import sys
import pandas as pd
import numpy as np
sys.path.append(os.path.join(proj_path, 'code'))
from utility import utility as util

geneanno_path = os.path.join(proj_path, 'data/curated/cancer.gene.anno.csv')
mutation_path = os.path.join(proj_path, 'data/UTSW_MW/All Mutations (CPRIT and COSMIC-84).csv')
out_path = os.path.join(proj_path, 'data/curated/lung/utsw.mw/lung_MutExpr_cancergene_loci.csv')

##############################    function    ##############################
def getMutationInfo(x):
    '''
    Get mutation position, mutation and source information
    '''
    try:
        source = re.findall(r'] \((.+)\)$', x)[0]
    except:
        sys.exit('Mutation source error: {}'.format(x))
    if source not in ['NCI-Navy', 'Gazdar']:
        try:
            pos = re.findall(r'[cn]\.([0-9-]+)', x.split(';')[0])[0]
        except:
            sys.exit('Mutation position error: {}'.format(x))
        try:
            mut = re.findall(r'[cn]\..*[0-9-]+(.*)', x.split(';')[0])[0]
        except:
            sys.exit('Mutation detail error: {}'.format(x))
        return pos, mut, source


################################    main    ###################################
mutation = pd.read_csv(mutation_path, low_memory=False)
mutation = mutation.drop(mutation.columns[[1,2]], axis=1)
mutation = pd.melt(mutation, id_vars=['Gene'], value_vars=mutation.columns[1:].tolist(), var_name='Cell', value_name='Mutation')
mutation = mutation.loc[mutation.Mutation.notnull(),:]
mutation.index = range(mutation.shape[0])

# parse mutation information
mut_list = []
for idx, row in mutation.iterrows():
    if idx % 100 == 1: print(idx)
    tmp_mut_tuple = [tmp_info for tmp_info in map(getMutationInfo, row['Mutation'].split(',')) if tmp_info is not None]
    if len(tmp_mut_tuple) > 0:
        tmp_pos, tmp_mut, tmp_source = zip(*tmp_mut_tuple)
        tmp_df = pd.DataFrame({'Gene':np.repeat(row['Gene'], len(tmp_pos)), 
                            'Cell':np.repeat(row['Cell'], len(tmp_pos)), 
                            'MutStart': tmp_pos,
                            'Mutation':tmp_mut,
                            'Source': tmp_source})
        mut_list.append(tmp_df)
mutation = pd.concat(mut_list)
print(mutation)