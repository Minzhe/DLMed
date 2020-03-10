########################################################################################
###                                curate.GeneSet.py                                 ###
########################################################################################

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import pandas as pd
import numpy as np
import pickle as pkl

##############################     function    #################################
def read_geneset(path, filtergenes=None, min_occ=5, min_perc=0.1, min_perc_small=0.5):
    '''
    Read GO term gene list
    '''
    geneset = dict()
    with open(path, 'r') as f:
        line = f.readline()
        while line:
            data = line.strip().split('\t')
            line = f.readline()
            name, genes = data[0], data[2:]
            if filtergenes is not None and min_perc is not None and min_occ is not None:
                # skip if all gene in the list is not contained in this gene set
                comm = set(genes) & set(filtergenes)
                if len(comm) >= min_occ:
                    if len(comm)/len(genes) >= min_perc:
                        geneset[name] = genes
                elif len(comm)/len(genes) > min_perc_small:
                    geneset[name] = genes
    return geneset

##############################     main    #################################
cancer_path = os.path.join(proj_dir, 'data/curated/cancer.gene.anno.csv')
c2_path = os.path.join(proj_dir, 'data/GeneSet/c2.cgp.v6.2.symbols.gmt')
c5_path = os.path.join(proj_dir, 'data/GeneSet/c5.all.v6.2.symbols.gmt')
out_path = os.path.join(proj_dir, 'data/curated/geneset.pkl')

cancergene = pd.read_csv(cancer_path, usecols=['Gene'], squeeze=True)
curated_geneset = read_geneset(c2_path, filtergenes=cancergene)
go_geneset = read_geneset(c5_path, filtergenes=cancergene)
all_geneset = {**curated_geneset, **go_geneset}
all_genes = list(set(gene for geneset in all_geneset.values() for gene in geneset))
cancergene = list(set(cancergene) & set(all_genes))
print(len(all_geneset), len(all_genes), len(cancergene))

with open(out_path, 'wb') as f:
    pkl.dump({'geneset':all_geneset, 'gene':sorted(all_genes), 'cancergene':sorted(cancergene)}, file=f)