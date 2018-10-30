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


def plot_density(title='', xlabel='', ylabel='', legend_title='', **data):
    '''
    Plot density curves.
    '''
    names = data.keys()
    colors = sns.color_palette(palette='hls', n_colors=len(names))
    f, ax = plt.subplots()
    for i, name in enumerate(names):
        seq = data[name]
        sns.distplot(seq, hist=False, kde=True, label=name, color=colors[i])
    ax.legend(title=legend_title)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return f, ax

def plot_multi_density(df, title='', xlabel='', xlim=None, ylim=None, xticks=[], yticks=[], oneplot=True, nrow=None, ncol=None):
    '''
    Plot density.
    '''
    if oneplot:
        f, ax = plt.subplots()
        for col in df:
            sns.distplot(df[col], hist=False, kde=True, color='steelblue', kde_kws={'alpha': 0.05}, ax=ax)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
    elif nrow is not None and ncol is not None:
        f, ax = plt.subplots(nrow, ncol)
        if df.shape[1] > nrow * ncol: raise ValueError('Unable to plot all columns because the number of rows and columns is not large enough.')
        num = 0
        for i in range(nrow):
            for j in range(ncol):
                sns.distplot(df.iloc[:,num], hist=False, kde=True, color='steelblue', ax=ax[i,j])
                ax[i,j].set_xlim(xlim)
                ax[i,j].set_ylim(ylim)
                ax[i,j].set_xticks(xticks)
                ax[i,j].set_yticks(yticks)
                ax[i,j].set_xlabel(df.columns[num], fontsize=8)
                num += 1
    else:
        raise ValueError('Incorrect number of rows and columns.')
    
    return f, ax