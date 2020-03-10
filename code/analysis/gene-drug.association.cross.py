####################################################################################################
###                                   gene-drug.association.py                                   ###
####################################################################################################
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
import os
import sys
sys.path.append(os.path.join(proj_dir, 'code'))
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#########################################    function   #########################################
def getPairRank(path):
    data = pd.read_csv(path, index_col=0)
    data.index.name = 'Gene'
    data = data.reset_index()
    data = pd.melt(data, id_vars=['Gene'], var_name='Drug', value_name='Effect')
    return list(data['Effect'].abs().sort_values(ascending=False).index)

def plotAgreement(s1, s2, path):
    n, i = len(s1), -(-len(s1) // 100)
    ratio = []
    interval = i
    while interval - i < n:
        x1, x2 = s1[:interval], s2[:interval]
        ratio.append(calOverlapRatio(x1, x2, n))
        interval += i
    ratio = pd.DataFrame(ratio, index=range(len(ratio)), columns=['Actual', 'Random'])
    # plot
    f, ax = plt.subplots(figsize=(8,8))
    plt.plot(ratio.index, ratio['Actual'])
    plt.plot(ratio.index, ratio['Random'])
    f.savefig(path)
    return ratio


def calOverlapRatio(x1, x2, total):
    actual = len(set(x1) & set(x2)) / len(set(x1) | set(x2))
    random = len(x1) / (2*total - len(x1))
    return actual, random


#########################################    main   #########################################
data_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.{}/{}.single_mut_simu.gene-drug.associ.median.gene.csv')
full_path = data_path.format('full_cnn', 'full_cnn')
gene_path = data_path.format('gene_cnn', 'gene_cnn')
img_path = os.path.join(proj_dir, 'result/simulation/single_mut_simulation_analysis.compare/full_gene.agreement.png')

full = getPairRank(full_path)
gene = getPairRank(gene_path)
ratio = plotAgreement(full, gene, img_path)
# print(ratio)