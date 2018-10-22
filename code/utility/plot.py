###################################################################################
###                                 plot.py                                     ###
###################################################################################

import seaborn as sns
import matplotlib.pyplot as plt
sns.set
sns.set_style('darkgrid')

def plotNegative(phenotable, path):
    '''
    Compare two methods for single gene knockout.
    '''
    same_idx = (phenotable.Gene1 == phenotable.Gene2) & (phenotable.Gene1 != 'negative')
    neg_idx = ((phenotable.Gene1 == 'negative') | (phenotable.Gene2 == 'negative')) & (phenotable.Gene1 != phenotable.Gene2)
    plt.clf()
    sns.distplot(phenotable.Score[same_idx], hist=False, kde=True, label='same')
    sns.distplot(phenotable.Score[neg_idx], hist=False, kde=True, label='negative')
    plt.legend(title = 'method')
    plt.title('Density plot of phenotype with single gene knock out')
    plt.xlabel('Phenotye score')
    plt.ylabel('Density')
    plt.savefig(path)


def plotDist(**scores):
    '''
    Plot density curves and negative line.
    '''
    names = scores.keys()
    colors = sns.color_palette(palette='hls', n_colors=len(names))
    plt.clf()
    for i, name in enumerate(names):
        score = scores[name]
        neg = score.loc[(score.iloc[:,0] == 'negative') & (score.iloc[:,1] == 'negative'),'score'].values
        sns.distplot(score.score.values.reshape(-1,1), hist=False, kde=True, label=name, color=colors[i])
        plt.axvline(neg, color=colors[i])
    plt.legend(title = 'method')
    plt.title('Density plot of scaled dual gene knock out phenotype score')
    plt.xlabel('Phenotye score')
    plt.ylabel('Density')
    plt.show()

def plot_density(df, xlim=None, ylim=None):
    '''
    Plot density.
    '''
    for col in df:
        sns.distplot(df[col], hist=False, kde=True, color='steelblue', kde_kws={'alpha': 0.05})
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel('')
