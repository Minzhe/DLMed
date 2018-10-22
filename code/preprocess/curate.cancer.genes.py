#################################################################################
###                         curate.cancer.genes.py                            ###
#################################################################################
# This script is to curate cancer related genes for further labelling.

import pandas as pd

########################   fucntion   ###########################
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
    


anno_path = "D:/projects/DLCell/data/Cosmic/Census_all.csv"
anno_out = "D:/projects/DLCell/data/Cosmic/cancer.gene.anno.csv"
uniprot_path = "D:/projects/DLCell/data/UniProt/uniprot.human.txt"

###########################   main   ###############################
### filter genes
gene_anno = pd.read_csv(anno_path, header=0)[['Gene Symbol', 'Tumour Types(Somatic)', 'Role in Cancer']].fillna('')
gene_anno.columns = ['Gene', 'Tissue', 'Role']
gene_anno.Role = gene_anno.apply(annotateGene, axis=1)
gene_anno = gene_anno[['Gene', 'Role']]

### prase uniprot
unip = pd.read_csv(uniprot_path, sep='\t')[['Entry', 'Entry name', 'Gene names']].fillna('')
unip.columns = ['Uniprot_AC', 'Uniprot_ID', 'Gene']
unip_dict = gene2Uniprot(unip, verbose=False)

### add uniprot protein name
gene_anno['Uniprot_AC'] = gene_anno.Gene.apply(lambda x: ';'.join(parseGene(unip_dict, x)[0]))
gene_anno['Uniprot_ID'] = gene_anno.Gene.apply(lambda x: ';'.join(parseGene(unip_dict, x)[1]))
print(gene_anno.head())

# gene_anno.to_csv(anno_out, index=None)
