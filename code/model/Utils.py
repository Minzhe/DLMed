##########################################################################
###                             Utils.py                               ###
##########################################################################
import numpy as np
import sys
import torch
import time
from itertools import combinations
from MGCNet import NCEAverage, NCESoftmaxLoss, NCECriterion

### ----------------------------  Dataset  ---------------------------- ###
class GeneticDataset(torch.utils.data.Dataset):
    __initialized = False
    def __init__(self, cells, data_dict):
        '''
        Args:
            cells: shared cells across different datasets
        '''
        self.cells = cells
        self.indexes = {c: i for i, c in enumerate(cells)}
        self.features = data_dict.keys()
        self.data_dict = data_dict
        self.__initialized = True

    def __len__(self):
        '''
        Denotes the number of samples
        '''
        return len(self.cells)
    
    def __getitem__(self, index):
        '''
        Generate one batch of data.
        Returns:
            idx: indexes of samples (long)
        '''
        # Generate indexes of the batch
        cell = self.cells[index]
        idx = self.indexes[cell]
        # Get data
        data = self.__data_generation(cell)
        return data, idx
    
    def __data_generation(self, indexes):
        '''
        Generates data containing batch_size samples.
        '''
        data = dict()
        for f in self.features:
            if indexes in self.data_dict[f].index:
                data[f] = torch.tensor(self.data_dict[f].loc[indexes, :].values)
            else:
                data[f] = self.__sample(self.data_dict[f], indexes)
        return data
    
    def __sample(self, df, indexes):
        mat = df.sample(n=1).values
        mat = np.random.permutation(mat.reshape(-1))
        return torch.tensor(mat)


class ContrastTrainer(object):
    '''
    Contrast learning
    '''
    def __init__(self, model, data_loader, contrast_param, optimizer, gradient_clip=10, device='cuda: 0'):
        self.model = model
        self.features = self.model.features
        self.data_loader = data_loader
        self.contrast_param = {'K': 100, 'T': 0.07, 'momentum': 0.9, 'use_softmax': False, 'device': device}
        self.contrast_param.update(contrast_param)
        contrast = lambda: NCEAverage(**self.contrast_param).to(device)
        self.contrasters = {(x, y): contrast() for (x, y) in combinations(self.features, 2)}
        criterion = lambda: NCESoftmaxLoss().to(device) if self.contrast_param['use_softmax'] else \
                            NCECriterion(self.contrast_param['output_size']).to(device)
        self.criterions = {(x, y, _): criterion() for (x, y) in self.contrasters.keys() for _ in [1,2]}
        self.optimizer = optimizer
        self.gradient_clip = gradient_clip
        self.device = device
        loss_hist = {(x, y, o, 'loss'): [] for (x, y, o) in self.criterions.keys()}
        prob_hist = {(x, y, o, 'prob'): [] for (x, y, o) in self.criterions.keys()}
        self.history = {**loss_hist, **prob_hist}
    
    def train(self, epochs, save_freq, model_path, hist_path):
        print('Start training ...')
        for epoch in range(1, epochs+1):
            self.train_one_epoch(epoch)
            if epoch % save_freq == 0:
                print('## Saving model to {}'.format(model_path.format(epoch)))
                state = {'model': self.model, 'optimizer': self.optimizer.state_dict(), 'epoch': epoch}
                for (f1, f2) in self.contrasters.keys():
                    state[('contrast', f1, f2)] = self.contrasters[(f1, f2)].state_dict()
                torch.save(state, model_path.format(epoch))
                print('## Saving history to {}'.format(hist_path))
                torch.save(self.history, hist_path)


    def train_one_epoch(self, epoch):
        # set train mode
        self.model.train()
        for cont in self.contrasters.keys():
            self.contrasters[cont].train()
        # set tracer
        batch_time = AverageMeter()
        total_loss = AverageMeter()
        tracer = {k: AverageMeter() for k in self.history.keys()}
        start_time = time.time()

        for idx, (X, index) in enumerate(self.data_loader):
            batch_size = X[list(X.keys())[0]].size(0)
            index = index.to(self.device)
            x = {key: val.float().to(self.device) for key, val in X.items()}

            # ---------------------- forward --------------------- #
            y = self.model(x)
            loss = 0
            for (f1, f2) in self.contrasters.keys():
                out1, out2 = self.contrasters[(f1, f2)](y[f1], y[f2], index)
                loss1 = self.criterions[(f1, f2, 1)](out1)
                loss2 = self.criterions[(f1, f2, 2)](out2)
                prob1 = out1[:,0].mean()
                prob2 = out2[:,0].mean()
                loss += loss1 + loss2
                # update tracer
                tracer[(f1, f2, 1, 'loss')].update(loss1.item(), batch_size)
                tracer[(f1, f2, 2, 'loss')].update(loss2.item(), batch_size)
                tracer[(f1, f2, 1, 'prob')].update(prob1.item(), batch_size)
                tracer[(f1, f2, 2, 'prob')].update(prob2.item(), batch_size)
            total_loss.update(loss.item(), batch_size)

            # ---------------------- backward --------------------- #
            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.gradient_clip)
            for key in self.contrasters.keys():
                torch.nn.utils.clip_grad_norm_(self.contrasters[key].parameters(), self.gradient_clip)
            self.optimizer.step()
            torch.cuda.synchronize()
            batch_time.update(time.time() - start_time)
        
        # -------------- print info ---------------- #
        print('Epoch: [{}]\t loss: {:.3f}\n'.format(epoch, total_loss.avg), end='')
        for (f1, f2) in combinations(self.features, 2):
            print('p_{}_{}: {:.3f},{:.3f}\t'.format(f1, f2, tracer[(f1, f2, 1, 'prob')].avg, tracer[(f1, f2, 2, 'prob')].avg), end='')
        print('\n------')
        sys.stdout.flush()

        # ----------------- debug ----------------- #
        if np.isnan(total_loss.avg):
            print(list(self.model.parameters()))
            print(X)
            raise ValueError('Nan detected.')

        # ------------- add to history ------------ #
        for (f1, f2, o, item) in self.history.keys():
            self.history[(f1, f2, o, item)].append(tracer[(f1, f2, o, item)].avg)




### ----------------------------  Help function  ---------------------------- ###
class AverageMeter(object):
    '''
    Computes and stores the average and current value
    '''
    def __init__(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count