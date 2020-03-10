###################################################################################
###                              curate.CTRP.py                                 ###
###################################################################################
# curate CTRP drug sensitivity dataset
# proj_dir = 'D:/projects/DLMed'
proj_dir = '/work/bioinformatics/s418336/projects/DLMed/'
import os
import pandas as pd

fit_path = os.path.join(proj_dir, 'data/CTRP/v20.data.curves_post_qc.txt')
comp_path = os.path.join(proj_dir, 'data/CTRP/v20.meta.per_compound.txt')
cell_path = os.path.join(proj_dir, 'data/CTRP/v20.meta.per_cell_line.txt')
exper_path = os.path.join(proj_dir, 'data/CTRP/v20.meta.per_experiment.txt')
out_path = os.path.join(proj_dir, 'data/curated/Lung/lung.CTRP.csv')

#######################    main    ############################
fit_data = pd.read_table(fit_path, sep='\t', usecols=['experiment_id', 'apparent_ec50_umol', 'master_cpd_id'])
comp_data = pd.read_table(comp_path, sep='\t', index_col=0, usecols=['master_cpd_id', 'cpd_name'], squeeze=True)
cell_data = pd.read_table(cell_path, sep='\t', index_col=0, usecols=['master_ccl_id', 'ccl_name', 'ccl_availability', 'ccle_primary_site'])
exper_data = pd.read_table(exper_path, sep='\t', index_col=0, usecols=['experiment_id', 'master_ccl_id'], squeeze=True).drop_duplicates()

# filter and map lung cancer cell line
cell_data.fillna('', inplace=True)
# cell_data = cell_data.loc[cell_data['ccle_primary_site'].str.upper().str.contains('LUNG'),:]

# map experiment to cell line
fit_data.experiment_id = fit_data.experiment_id.apply(lambda x: exper_data.loc[x] if x in exper_data.index else None)
fit_data = fit_data.loc[fit_data.experiment_id.notnull(),:]
fit_data.experiment_id = fit_data.experiment_id.apply(lambda x: cell_data.loc[x,'ccl_name'] if x in cell_data.index else None)
fit_data = fit_data.loc[fit_data.experiment_id.notnull(),:]

# map drug id to name
fit_data.master_cpd_id = fit_data.master_cpd_id.apply(lambda x: comp_data.loc[x])

# clean name
fit_data.columns = ['Cell', 'EC50', 'Drug']
fit_data = fit_data[['Cell', 'Drug', 'EC50']]
fit_data.sort_values(by=['Cell', 'Drug'], inplace=True)
print(fit_data)
exit()

fit_data.to_csv(out_path, index=None)