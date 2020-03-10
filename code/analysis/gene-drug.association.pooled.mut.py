###############################################################################################################
###                                   gene-drug.association.pooled.mut.py                                   ###
###############################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
import pandas as pd
import numpy as np
import pickle as pkl

###################################    function   ######################################
def poolGeneDrugAssoci(mut, res):
    '''
    Calculate the gene drug interaction by simple two group comparison
    '''
    # count mutation allele frequency
    mut_freq = mut.sum()
    mut_freq = mut_freq[mut_freq != 0].apply(int)
    # calculate the loci drug response effect
    eff = pd.DataFrame(np.nan, index=mut_freq.index, columns=res.columns)
    eff.columns.name = None
    # get group index for each mutation loci
    for loci in eff.index:
        print(loci)
        no_mut_res = res.loc[mut.index[mut[loci] == 0],:].median()
        mut_res = res.loc[mut.index[mut[loci] == 1],:].median()
        eff.loc[loci,:] = mut_res - no_mut_res
    eff.insert(loc=0, column='Freq', value=mut_freq)
    return eff
    

#####################################    main   #########################################
model = 'gene_cnn'
mut_loci_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mut_cancergene_loci_table.csv')
mut_gene_path = os.path.join(proj_dir, 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mut_cancergene_table.csv')
res_path = os.path.join(proj_dir, 'result/prediction/{}/results.csv'.format(model))
actual_loci_path = os.path.join(proj_dir, 'result/simulation/single_mut_pooled_analysis/single_mut_pooled.gene-drug.association.actual.loci.csv')
actual_gene_path = os.path.join(proj_dir, 'result/simulation/single_mut_pooled_analysis/single_mut_pooled.gene-drug.association.actual.gene.csv')
pred_path = os.path.join(proj_dir, 'result/simulation/single_mut_pooled_analysis/single_mut_pooled.gene-drug.association.pred.{}.csv'.format(model))

mut = pd.read_csv(mut_gene_path, index_col=0)
res = pd.read_csv(res_path)

res_actual = res.drop(['ic50_hat'], axis=1).pivot(index='cellline_name', columns='drug_name', values='ic50')
actual = poolGeneDrugAssoci(mut, res_actual)
actual.to_csv(actual_gene_path)

# res_pred = res.drop(['ic50'], axis=1).pivot(index='cellline_name', columns='drug_name', values='ic50_hat')
# pred = poolGeneDrugAssoci(mut, res_pred)
# pred.to_csv(pred_path)
