#########################################################################################
###                                curate.mutation.py                                 ###
#########################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os
import glob
import pandas as pd
import numpy as np

luad_path = os.path.join(proj_dir, 'data/TCGA/LUAD.Mutation_Packager_Oncotated_Calls')
lusc_path = os.path.join(proj_dir, 'data/TCGA/LUSC.Mutation_Packager_Oncotated_Calls')
geneanno_path = os.path.join(proj_dir, 'data/COSMIC/Census_all.csv')
out_path = os.path.join(proj_dir, 'data/curated/Lung/tcga/Lung.mutation.csv')

#################################    mutation   #################################
gene_list = list(pd.read_csv(geneanno_path, usecols=[0], squeeze=True))

mutation = list()
for sub_type in ['LUAD', 'LUSC']:
    if sub_type == 'LUAD':
        files = glob.glob(os.path.join(luad_path, '*.maf.txt'))
    elif sub_type == 'LUSC':
        files = glob.glob(os.path.join(lusc_path, '*.maf.txt'))
    for file in files:
        print('Reading {}'.format(file))
        mut = pd.read_table(file, sep='\t', skiprows=3, usecols=['Hugo_Symbol', 'Chromosome', 'Start_position', 'Variant_Classification', 'Tumor_Sample_Barcode'])
        mut = mut.loc[mut['Hugo_Symbol'].isin(gene_list),:]
        mut = mut.loc[mut['Variant_Classification'] != 'Silent',:]
        mut.columns = ['Gene', 'Chr', 'Mut_Start', 'Varient_Type', 'Patient']
        mut = mut[['Patient', 'Gene', 'Chr', 'Mut_Start', 'Varient_Type']]
        mut.sort_values(by=['Gene', 'Chr', 'Mut_Start'], inplace=True)
        mut['Type'] = sub_type
        mutation.append(mut)
mutation = pd.concat(mutation, axis=0)
mutation.to_csv(out_path, index=None)
