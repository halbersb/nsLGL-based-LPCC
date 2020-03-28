import numpy as np
import pandas as pd
import os

from regain.covariance import LatentTimeGraphicalLasso
from regain.datasets import make_dataset
from regain.utils import error_norm_time
from numpy import genfromtxt

np.random.seed(42)
pathX = 'databases\\G1_forLTGL\\X\\'
pathY = 'databases\\G1_forLTGL\\Y\\'
pathSave = 'databases\\After_LTGL\\'
for i in range (2):
    folders = os.listdir(pathX + 'G' + str(i+1))
    for item in folders:
        for j in range(5):
            print(pathX + 'G' + str(i+1) + '\\' + item + '\\' + item + '_' + str(j+1) + '.csv')
            X = genfromtxt(pathX + 'G' + str(i+1) + '\\' + item + '\\' + item + '_' + str(j+1) + '.csv', delimiter=',')
            Y = genfromtxt(pathY + 'G' + str(i+1) + '\\' + item + '\\' + item + '_' + str(j+1) + '.csv', delimiter=',')
            Y = Y.astype(int)
            mdl = LatentTimeGraphicalLasso(max_iter=50).fit(X, Y)
            results = np.empty([mdl.precision_.shape[1], mdl.precision_.shape[2]])
            for p in range(mdl.precision_.shape[0]):
                results = np.append(results, mdl.precision_[p, :, :], axis=0)
            cur_pathSave=pathSave + 'G' + str(i+1) + '\\' + item
            if not os.path.exists(cur_pathSave):
                os.makedirs(cur_pathSave)
            pd.DataFrame(results).to_csv(cur_pathSave + '\\' + item  + '_' + str(j) + '.csv')