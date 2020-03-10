####################################################################################
###                             train_MGCNet.py                                  ###
####################################################################################
import pandas as pd
import numpy as np
import os
import sys
import torch
import matplotlib.pyplot as plt
proj_dir = '/work/bioinformatics/s418336/projects/DLMed'
os.chdir(proj_dir)
sys.path.append('code/model')
from Utils import GeneticDataset, AverageMeter, ContrastTrainer
from MGCNet import MGCNet, NCEAverage, NCECriterion, NCESoftmaxLoss


##########################     functions    ##############################
def plot_hist(hist, fig_path):
    '''
    Plot training log
    '''
    fig, ax = plt.subplots(2, 1, figsize=(15,15))
    for key, val in hist.items():
        if 'loss' in key:
            ax[0].plot(np.arange(len(val)), val, label=key)
        if 'prob' in key:
            ax[1].plot(np.arange(len(val)), val, label=key)
    ax[0].set_xlabel('Epochs')
    ax[1].set_xlabel('Epochs')
    ax[0].set_ylabel('Loss')
    ax[1].set_ylabel('Prob')
    ax[0].legend()
    ax[1].legend()
    fig.savefig(fig_path)

def shuffle_df(df):
    '''
    Shuffle data frame
    '''
    shape = df.values.shape
    mat = np.random.permutation(df.values.reshape(-1))
    return pd.DataFrame(mat.reshape(shape), index=df.index, columns=df.columns)


###########################       main       ##############################

### ---------  1. load data  ---------- ###
mut_path = 'data/curated/Lung/merge_final_version/ccle_utsw.lung_Mut_cancergene_table.csv'
expr_path = 'data/curated/Lung/merge_final_version/ccle_utsw.lung_RNAseq_cancergene.csv'
cnv_path = 'data/curated/Lung/merge_final_version/ccle_utsw.lung_CNV_cancergene_table.csv'
drug_path = 'data/curated/Lung/merge_final_version/ccle_utsw.lung_drug.split.csv'

genetic_model_path = 'code/predict/model_repo/MGCNet.geneitc_embedding.shuffle_batch.epoch_{}.pt'
hist_path = 'code/predict/model_repo/MGCNet.hist.shuffle_batch.pt'
fig_path = 'code/predict/model_repo/MGCNet.log.shuffle_batch.png'

mut = pd.read_csv(mut_path, index_col=0)
expr = pd.read_csv(expr_path, index_col=0)
cnv = pd.read_csv(cnv_path, index_col=0)
# mat = shuffle_df(mut)
# expr = shuffle_df(expr)
# cnv = shuffle_df(cnv)
# drug = pd.read_csv(drug_path)
# drug = drug.loc[~drug['Dropped_previously'],:]

cells = sorted(list(set(expr.index)))
print('Unique cell lines: {}'.format(len(cells)))

### ---------- 2. data loader ----------- ###
train_data = GeneticDataset(cells, {'mut': mut, 'expr': expr, 'cnv': cnv})
train_loader = torch.utils.data.DataLoader(train_data, batch_size=8, shuffle=True, num_workers=1)

### ------------ 3. set MGCNet model ---------------- ###
device = 'cuda: 0'
n_data = len(cells)
n_out_features = 64

# NCE parameters
# Increasing nce_m improves stability
# !!! Adding batch normalization layer improves stability 
contrast_param = {'input_size': n_out_features, 'output_size': n_data, 'K': 100, 'T': 0.07, 'momentum': 0.9, 'use_softmax': False}

# Set model
genetic_model = MGCNet(layer_sizes={'mut': [mut.shape[1], 128, 128], 
                                    'expr': [expr.shape[1], 128, 128],
                                    'cnv': [cnv.shape[1], 128, 128]}, 
                       n_out_features=n_out_features).to(device)

# Set optimizer
gradient_clip = 5
# optimizer = torch.optim.SGD(genetic_model.parameters(), lr=1e-3, momentum=0.9, weight_decay=1e-4)
optimizer = torch.optim.Adam(genetic_model.parameters(), lr=1e-4)

### ---------------- 4. train MGCNet model ----------------- ###
# contrast_train = ContrastTrainer(genetic_model, train_loader, contrast_param, optimizer, gradient_clip, device)
# contrast_train.train(500, save_freq=20, model_path=genetic_model_path, hist_path=hist_path)

### ---------------- 4. plot history ----------------- ###
hist = torch.load(hist_path)
print('Plotting training history ...')
plot_hist(hist, fig_path)

### ----------------- 5. get embedding ---------------- ###
# checkpoint = torch.load(genetic_model_path.format(2800), map_location=device)
# genetic_model.load_state_dict(checkpoint['model'])
# for X, idx in train_loader:
#     for _ in X.keys():
#         X[_] = X[_].float().to(device)
#     y = genetic_model(X)
#     print(y)
#     exit()