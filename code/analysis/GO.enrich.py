########################################################################################
###                                   GO.enrich.py                                   ###
########################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import glob
import pandas as pd
import numpy as np
import re
from scipy.stats import percentileofscore
import pickle as pkl
import multiprocessing as mp
from utility import utility as util
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from utility import plot as p

##############################    main   #################################
def read_simu(paths, loci2gene, cell_index, drug_index, perc_cutoff, out_dir):
    '''
    Load simulation result
    '''
    enriched_genes = dict()
    for path in list(paths):
        print('Loading {}'.format(path))
        with open(path, 'rb') as f:
            data = pkl.load(f)
            perb_ic50, base_ic50 = data['res'][0], np.array(data['base'][0]['y_hat_me_cnn'])
            cell, drugs = cell_index[data['cellline_id']], drug_index.values()
        diff_logic50 = perb_ic50-base_ic50
        resis_idx = diff_logic50 > np.log(1+perc_cutoff)
        sens_idx = diff_logic50 < np.log(1-perc_cutoff)
        diff_logic50[resis_idx | sens_idx] = 1
        diff_logic50[~(resis_idx | sens_idx)] = 0
        diff_logic50 = pd.DataFrame(diff_logic50, index=loci2gene, columns=drugs)
        diff_logic50 = diff_logic50.groupby(by=diff_logic50.index).any()
        enriched_genes.update({(cell,drug): list(diff_logic50.index.values[diff_logic50[drug] == 1]) for drug in diff_logic50.columns})
    out_path = os.path.join(out_dir, 'single_eff_mut.cutoff_{}.pkl'.format(perc_cutoff))
    with open(out_path, 'wb') as f:
        pkl.dump(enriched_genes, file=f)
    return enriched_genes


def enrich_test(enriched_genes, genesets, poolgenes):
    '''
    Do enrichment analysis for each drug in each cell line, and then calculate the percentage of cell line is enriched in the pathway respect to the drug.
    '''
    cells = np.unique(list(zip(*enriched_genes.keys()))[0])
    drugs = np.unique(list(zip(*enriched_genes.keys()))[1])
    enrich_sets = pd.DataFrame(0, index=genesets.keys(), columns=drugs)
    for drug in drugs:
        for go, geneset in genesets.items():
            # multiprocess training
            manager = mp.Manager()
            p_val = manager.dict()
            jobs = []
            for cell in cells:
                p = mp.Process(target=util.hyper_test, args=(enriched_genes[(cell, drug)], geneset, poolgenes, False, p_val, (go,drug,cell),))
                jobs.append(p)
                p.start()
            for proc in jobs:
                proc.join()
            p_val = np.array(p_val.values())
            enrich_sets.loc[go,drug] = round(sum(p_val<=0.05)/len(p_val), 4)
            print('Analyzing drug {} for pathway {}: {}'.format(drug, go, enrich_sets.loc[go,drug]))
    # filter 0 row and columns
    enrich_sets = enrich_sets.loc[enrich_sets.sum(axis=0)!=0, enrich_sets.sum(axis=1)!=0]
    return enrich_sets


def enrich_test_pool(enriched_genes, genesets, poolgenes):
    '''
    Do enrichment analysis for each drug in each pathway, multiple cell line are pooled.
    '''
    cells = np.unique(list(zip(*enriched_genes.keys()))[0])
    drugs = np.unique(list(zip(*enriched_genes.keys()))[1])
    enrich_sets = pd.DataFrame(np.nan, index=genesets.keys(), columns=drugs)
    # multiprocess training
    manager = mp.Manager()
    p_val = manager.dict()
    jobs = list()
    for drug in drugs:
        print(drug)
        tmp_genes = list(set(gene for cell in cells for gene in enriched_genes[(cell, drug)]))
        p = mp.Process(target=util.hyper_tests, args=(tmp_genes, genesets, poolgenes, False, p_val, drug))
        jobs.append(p)
        p.start()
    for proc in jobs:
        proc.join()
    
    p_val = dict(p_val)
    enrich_sets = pd.DataFrame(p_val)
    # filter 0 row and columns
    enrich_sets = enrich_sets.loc[enrich_sets.sum(axis=1)!=0, enrich_sets.sum(axis=0)!=0]
    return enrich_sets


def clean_name(x):
    x = x.replace('SIGNALING_', '')
    x = x.replace('APOPTOTIC_PATHWAY', 'APOPTOSIS')
    x = x.replace('POSITIVE_REGULATION', 'POS_REG')
    x = x.replace('NEGATIVE_REGULATION', 'NEG_REG')
    x = re.sub(r'_INVOLVED.*', '', x)
    x = re.sub(r'_VIA.*', '', x)
    x = re.sub(r'_BY_.*', '', x)
    x = re.sub(r'_INDUCED.*', '', x)
    x = re.sub(r'_IN_.*', '', x)
    # x = re.sub(r'_RESPONSE.*', '', x)
    return x

##############################    main   #################################
geneset_path = os.path.join(proj_dir, 'data/curated/geneset.pkl')
simu_path = glob.glob(os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation/ic50_cnn_model_cell_*_drug_*.pkl'))
loci_anno_path = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation/ic50_cnn_model_alleles.pkl')
simu_out_dir = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation')
simu_enrich_pool_path = os.path.join(proj_dir, 'result/simulation/random_single_mut_simulation.enrichgeneset_pool.csv')
heat_out_path = os.path.join(proj_dir, 'result/simulation/fig/enrichgeneset_pool.cancer_gene.heatmap.png')
geneset_names = ['Notch', 'Wnt', 'HIF', 'GO_APOPTOTIC', 'EXTRINSIC_APOPTOTIC', 'INTRINSIC_APOPTOTIC', 'TGFB', 'P53', 'AKT',
                 'cell_cycle', 'cell_death', 'JAK_STAT', 'MAPK', 'PPAR', 'ERBB',
                 'KRAS', 'GO_CELL_PROLIFERATION', 'cell_growth', 'chemo', 'drug', 
                 'lung_cancer']

# geneset
with open(geneset_path, 'rb') as f:
    data = pkl.load(f)
    all_genesets, all_genes, cancergene = data['geneset'], data['gene'], data['cancergene']
# location annotation
with open(loci_anno_path, 'rb') as f:
    data = pkl.load(f)
    loci, gene_index, drug_index, cell_index = data['alleles'], data['gene_index'], data['drug_index'], data['cell_index']
    loci2gene = [gene_index[l].split('_')[0] for l, _ in loci]
# simulated enriched gene
# eff_genes = read_simu(simu_path, loci2gene=loci2gene, cell_index=cell_index, drug_index=drug_index, perc_cutoff=0.15, out_dir=simu_out_dir)
with open(os.path.join(simu_out_dir, 'single_eff_mut.cutoff_0.15.pkl'), 'rb') as f:
    enriched_genes = pkl.load(f)

# calculate pooled p-value
# enrich_sets = enrich_test_pool(enriched_genes, genesets=all_genesets, poolgenes=cancergene)
enrich_sets = pd.read_csv(simu_enrich_pool_path, index_col=0)
# f, ax = plt.subplots(figsize=(20,18))
geneset_names = [set_.upper() for set_ in geneset_names]
enrich_sets = enrich_sets.loc[enrich_sets.index.to_series().apply(lambda x: any([name in x for name in geneset_names])),:]
enrich_sets.index = enrich_sets.index.to_series().apply(clean_name)
f, ax = plt.subplots(figsize=(40,30))
p.plotHeatmap(enrich_sets, cluster=True, ax=ax)
f.savefig(heat_out_path, transparent=True)
