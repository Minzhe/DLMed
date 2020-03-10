##############################################################################
###                             MGCNet.py                                  ###
##############################################################################
# Refer to CMC: https://github.com/HobbitLong/CMC/blob/master/train_CMC.py
import math
import torch
from torch import nn
import torch.nn.functional as F

####################               network                #######################
class MGCNet(torch.nn.Module):
    '''
    Network for contrastive coding of multiple genetic information
    '''
    def __init__(self, layer_sizes, n_out_features):
        '''
        Args:
            layer_sizes: e.g., {"mutation": [10000, 100, 100]} to define the structures, the keys should exist in keys for DataLoader
            n_out_features: number of output features
        '''
        super(MGCNet, self).__init__()
        self.features = list(layer_sizes.keys())
        self.out_shape = n_out_features
        for f in self.features:
            net = []
            in_shape = layer_sizes[f][0]
            for i, n_neurons in enumerate(layer_sizes[f][1:]):
                net.append(torch.nn.Linear(in_shape, n_neurons))
                net.append(torch.nn.BatchNorm1d(n_neurons))
                net.append(torch.nn.ReLU6())
                in_shape = n_neurons
            net.append(torch.nn.Linear(in_shape, self.out_shape))
            setattr(self, f, torch.nn.ModuleList(net))
        
    def forward(self, X):
        '''
        Args:
            x: a dict
        '''
        out = dict()
        for f in self.features:
            x_ = X[f]
            for layer in getattr(self, f):
                x_ = layer(x_)
            out[f] = x_
        return out


# class neural_net(torch.nn.Module):
#     '''
#     Network to predict drug IC50
#     '''
#     def __init__(self, input_shapes, layer_sizes):
#         super().__init__()
#         net = []
#         in_shape = sum(input_shapes)
#         for i, n_neurons in enumerate(layer_sizes):
#             net.append(torch.nn.Linear(in_shape, n_neurons))
#             net.append(torch.nn.BatchNorm1d(n_neurons))
#             net.append(torch.nn.ReLU6())
#             in_shape = n_neurons
#         net.append(torch.nn.Linear(in_shape, 1))
#         self.net = torch.nn.ModuleList(net)
    
#     def forward(self, X):



####################              loss                #######################
class NCEAverage(nn.Module):
    def __init__(self, input_size, output_size, K, T=0.07, momentum=0.5, use_softmax=False, device='cuda: 0'):
        '''
        Args:
            input_size: n_features
            output_size: n_samples
            K: number of negatives to constrast for each positive
            T: temperature that modulates the distribution
        '''
        super(NCEAverage, self).__init__()
        self.output_size = output_size
        self.unigrams = torch.ones(self.output_size)
        self.multinomial = AliasMethod(self.unigrams)
        self.multinomial.to(device)
        self.K = K
        self.use_softmax = use_softmax

        self.register_buffer('params', torch.tensor([K, T, -1, -1, momentum]))
        stdv = 1. / math.sqrt(input_size / 3)
        self.register_buffer('memory_x1', torch.rand(output_size, input_size).mul_(2 * stdv).add_(-stdv))
        self.register_buffer('memory_x2', torch.rand(output_size, input_size).mul_(2 * stdv).add_(-stdv))

    def forward(self, x1, x2, index, idx=None):
        '''
        Args:
            x1: out_features for x1
            x2: out_features for x2
            index: torch.long for data ids
        '''
        K = int(self.params[0].item())
        T = self.params[1].item()

        momentum = self.params[4].item()
        batch_size = x1.size(0)
        input_size = self.memory_x1.size(1)

        # score computation
        if idx is None:
            idx = self.multinomial.draw(batch_size * (self.K + 1)).view(batch_size, -1)
            idx.select(1, 0).copy_(index.data)
        # sample
        weight_x2 = torch.index_select(self.memory_x2, 0, idx.view(-1)).detach()
        weight_x2 = weight_x2.view(batch_size, K + 1, input_size)
        out_x1 = torch.bmm(weight_x2, x1.view(batch_size, input_size, 1))
        # sample
        weight_x1 = torch.index_select(self.memory_x1, 0, idx.view(-1)).detach()
        weight_x1 = weight_x1.view(batch_size, K + 1, input_size)
        out_x2 = torch.bmm(weight_x1, x2.view(batch_size, input_size, 1))
        # Batchwise matrix multiplication
        # weight_x1: [batch_size, K + 1, n_out_features]
        # x2:        [batch_size, n_out_features]
        # out_x2:    [batch_size, K + 1, 1]
        
        if self.use_softmax:
            out_x1 = torch.div(out_x1, T)
            out_x2 = torch.div(out_x2, T)
            out_x1 = out_x1.contiguous()
            out_x2 = out_x2.contiguous()
        else:
            out_x1_e = torch.exp(out_x1 - torch.max(out_x1, dim=1, keepdim=True)[0])
            out_x1_s = torch.sum(out_x1_e, dim=1, keepdim=True)
            out_x1 = torch.div(out_x1_e, out_x1_s)
            out_x2_e = torch.exp(out_x2 - torch.max(out_x1, dim=1, keepdim=True)[0])
            out_x2_s = torch.sum(out_x2_e, dim=1, keepdim=True)
            out_x2 = torch.div(out_x2_e, out_x2_s)

        # # update memory
        with torch.no_grad():
            x1_pos = torch.index_select(self.memory_x1, 0, index.view(-1))
            x1_pos.mul_(momentum)
            x1_pos.add_(torch.mul(x1, 1 - momentum))
            x1_norm = x1_pos.pow(2).sum(1, keepdim=True).pow(0.5)
            updated_x1 = x1_pos.div(x1_norm)
            self.memory_x1.index_copy_(0, index, updated_x1)
            x2_pos = torch.index_select(self.memory_x2, 0, index.view(-1))
            x2_pos.mul_(momentum)
            x2_pos.add_(torch.mul(x2, 1 - momentum))
            x2_norm = x2_pos.pow(2).sum(1, keepdim=True).pow(0.5)
            updated_x2 = x2_pos.div(x2_norm)
            self.memory_x2.index_copy_(0, index, updated_x2)

        return out_x1, out_x2

eps = 1e-7
class NCECriterion(nn.Module):
    """
    Eq. (12): L_{NCE}
    """
    def __init__(self, n_data):
        super(NCECriterion, self).__init__()
        self.n_data = n_data

    def forward(self, x):
        bsz = x.shape[0]
        m = x.size(1) - 1

        # noise distribution
        Pn = 1 / float(self.n_data)

        # loss for positive pair
        P_pos = x.select(1, 0)
        log_D1 = torch.div(P_pos, P_pos.add(m * Pn + eps)).log_()

        # loss for K negative pair
        P_neg = x.narrow(1, 1, m)
        log_D0 = torch.div(P_neg.clone().fill_(m * Pn), P_neg.add(m * Pn + eps)).log_()

        loss = - (log_D1.sum(0) + log_D0.view(-1, 1).sum(0)) / bsz

        return loss


class NCESoftmaxLoss(nn.Module):
    """Softmax cross-entropy loss (a.k.a., info-NCE loss in CPC paper)"""
    def __init__(self):
        super(NCESoftmaxLoss, self).__init__()
        self.criterion = nn.CrossEntropyLoss()

    def forward(self, x):
        bsz = x.shape[0]
        x = x.squeeze()
        label = torch.zeros([bsz]).cuda().long()
        loss = self.criterion(x, label)        
        print(loss)        
        return loss


####################           function           ####################
class AliasMethod(object):
    """
    From: https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
    """
    def __init__(self, probs):

        if probs.sum() > 1:
            probs.div_(probs.sum())
        K = len(probs)
        self.prob = torch.zeros(K)
        self.alias = torch.LongTensor([0]*K)

        # Sort the data into the outcomes with probabilities
        # that are larger and smaller than 1/K.
        smaller = []
        larger = []
        for kk, prob in enumerate(probs):
            self.prob[kk] = K*prob
            if self.prob[kk] < 1.0:
                smaller.append(kk)
            else:
                larger.append(kk)

        # Loop though and create little binary mixtures that
        # appropriately allocate the larger outcomes over the
        # overall uniform mixture.
        while len(smaller) > 0 and len(larger) > 0:
            small = smaller.pop()
            large = larger.pop()

            self.alias[small] = large
            self.prob[large] = (self.prob[large] - 1.0) + self.prob[small]

            if self.prob[large] < 1.0:
                smaller.append(large)
            else:
                larger.append(large)

        for last_one in smaller+larger:
            self.prob[last_one] = 1

    def cuda(self):
        self.prob = self.prob.cuda()
        self.alias = self.alias.cuda()
        
    def to(self, device):
        self.prob = self.prob.to(device)
        self.alias = self.alias.to(device)

    def draw(self, N):
        """
        Draw N samples from multinomial
        :param N: number of samples
        :return: samples
        """
        K = self.alias.size(0)

        kk = torch.zeros(N, dtype=torch.long, device=self.prob.device).random_(0, K)
        prob = self.prob.index_select(0, kk)
        alias = self.alias.index_select(0, kk)
        # b is whether a random number is greater than q
        b = torch.bernoulli(prob)
        oq = kk.mul(b.long())
        oj = alias.mul((1-b).long())

        return oq + oj
