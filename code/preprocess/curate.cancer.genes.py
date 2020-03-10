#################################################################################
###                         curate.cancer.genes.py                            ###
#################################################################################
# This script is to curate cancer related genes for further labelling.

import os
import numpy as np
import pandas as pd
from functools import reduce

proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
cancer_gene_path = os.path.join(proj_dir, 'data/GeneSet/Census_all.csv')
oncogenic_path = os.path.join(proj_dir, 'data/GeneSet/c6.all.v6.2.symbols.gmt')
cgn_path = os.path.join(proj_dir, 'data/GeneSet/c4.cgn.v6.2.symbols.gmt')
pertub_path = os.path.join(proj_dir, 'data/GeneSet/c2.cgp.v6.2.symbols.gmt')
geneset_out = os.path.join(proj_dir, 'data/curated/cancer.gene.csv')
uniprot_path = os.path.join(proj_dir, 'data/UniProt/uniprot.human.txt')

########################   fucntion   ###########################
def read_gsea(file):
    '''
    Read GESA gene set .gmt file.
    '''
    with open(file) as f:
        lines = f.readlines()
    
    name, url, genes = list(), list(), list()
    for line in lines:
        items = line.strip().split('\t')
        name.append(items[0])
        url.append(items[1])
        genes.append(';'.join(items[2:]))

    return pd.DataFrame({'Name': name, 'Url': url, 'Genes': genes})     


def annotateGene(row):
    '''
    Annotate genes based on census description of roles in cancer.
    Annotation of some gene with both oncogenic and tumor suppressive role:
    ---------
    ARID1A: TSG in breast cancer
    EP300: TSG in breast cancer
    ESR1: TSG and oncogene
    ETV6: oncogene breast cancer
    GATA3: TSG and oncogene
    IRS4: oncogene
    MALAT1: oncogene
    MAP2K4: oncogene
    MAP3K1: oncogene
    MAP3K13: oncogene
    NOTCH1: TSG and oncogene
    RAD21: oncogene in breast cancer
    TBX3: oncogene
    TP53: TSG
    --------
    '''
    gene, tissue, role = row['Gene'], row['Tissue'], row['Role'].replace('fusion', 'oncogene')

    if 'TSG' in role and 'oncogene' in role:
        if 'breast' in tissue or 'lung' in tissue:
            print('Warning: {}\tboth TSG and oncogene.'.format(gene))
        return 'unclear'
    elif 'TSG' in role:
        return 'TSG'
    elif 'oncogene' in role:
        return 'oncogene'
    elif role == '':
        return 'unknown'
    else:
        raise ValueError('Unparsed description: {}'.format(role))


def gene2Uniprot(anno, verbose=True):
    '''
    Create dictionary to store gene to uniprot id information.
    '''
    unip_dict = dict()
    for idx, row in anno.iterrows():
        genes = row['Gene'].split(' ')
        for gene in genes:
            if gene in unip_dict.keys():
                if verbose:
                    print('Warning: Gene already exist in dictionary.\nOld: {}, {}.\nNow: {}, {}.'.format(gene, unip_dict[gene], gene, (row['Uniprot_AC'], row['Uniprot_ID'])))
                unip_dict[gene] = unip_dict[gene][0] + [row['Uniprot_AC']], unip_dict[gene][1] + [row['Uniprot_ID']]
            else:
                unip_dict[gene] = [row['Uniprot_AC']], [row['Uniprot_ID']]
    
    return unip_dict

def editDict(genedict):
    '''
    '''
    genedict['IGK'] = genedict['IGK@']
    genedict['DUXL4']
    del genedict['IGK@']
    

###########################   main   ###############################
### census genes
gene_anno = pd.read_csv(cancer_gene_path, header=0)[['Gene Symbol', 'Tumour Types(Somatic)', 'Role in Cancer']].fillna('')
gene_anno.columns = ['Gene', 'Tissue', 'Role']
gene_anno.Role = gene_anno.apply(annotateGene, axis=1)
gene_anno = gene_anno[['Gene', 'Role']]
cancer_gene = list(gene_anno['Gene'])

### GSEA oncogenic gene set
oncogenic = read_gsea(oncogenic_path)
oncogenic_gene = reduce(lambda a,b: a+b, oncogenic.Genes.apply(lambda x: x.split(';')))
oncogenic_gene = np.unique(pd.Series(oncogenic_gene))
print(oncogenic.shape[0])
print(oncogenic_gene)
print(len(oncogenic_gene))
print(len(set(cancer_gene) | set(oncogenic_gene)))

### GSEA cancer gene neighborhoods
cgn = read_gsea(cgn_path)
cgn_gene = reduce(lambda a,b: a+b, cgn.Genes.apply(lambda x: x.split(';')))
cgn_gene = np.unique(pd.Series(cgn_gene))
print(cgn.shape[0])
print(cgn_gene)
print(len(cgn_gene))
print(len(set(cancer_gene) | set(cgn_gene)))

### GSEA chemical pertubation
pertub = read_gsea(pertub_path)
pertub_gene = reduce(lambda a,b: a+b, pertub.Genes.apply(lambda x: x.split(';')))
pertub_gene = np.unique(pertub_gene)
print(pertub.shape[0])
print(pertub_gene)
print(len(pertub_gene))
print(len(set(cancer_gene) | set(pertub_gene)))

