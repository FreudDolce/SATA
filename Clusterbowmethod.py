# ===============================
# Used for patient cluster
# 20200920
# Ji Hongchen
# ===============================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
import os
import math
from scipy import stats
from scipy.spatial.distance import cdist
from AMBscattersurvivalforgene import SurvivalAnalysis, SurvivalScatter

MEAND = []


CMAP = 'Greens'
N_CLUSTER = 7
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_patient_20200617/'
OUTPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/cluster/' + \
    str(N_CLUSTER) + '/'
if os.path.exists(OUTPATH) == False:
    os.mkdir(OUTPATH)
CLUSTER_COL = []
CLASS_LIST = []
for i in range(10):
    CLUSTER_COL.append('RateClass_' + str(i + 1))
for j in range(N_CLUSTER):
    CLASS_LIST.append('Class_' + str(j + 1))


def load_data(dataframe):
    dataframe['to_death_time'][dataframe['to_death_time']. isin(
        ['data_loss', 'data_empty'])] = 0
    dataframe['to_last_followup'][dataframe['to_last_followup'].isin(['data_loss',
                                                                      'data_empty'])] = 0
    for item in ['to_death_time', 'to_last_followup']:
        dataframe[item] = pd.to_numeric(dataframe[item], errors='coerce')
    dataframe['to_last_known_alive'] = dataframe['to_death_time'] + \
        dataframe['to_last_followup']
    dataframe['patient_status'][dataframe['patient_status'] == 'Alive'] = 0
    dataframe['patient_status'][dataframe['patient_status'] == 'Dead'] = 1
    dataframe.drop(dataframe[dataframe['patient_status'].isin(['data_loss',
                                                               'data_empty'])].index, inplace=True)


def KmeansClusterTrain(trainframe, n_k):
    traindata = np.array(trainframe[CLUSTER_COL])
    kmc = KMeans(n_clusters=n_k, max_iter=3000,
                 n_init=100, init='k-means++')
    kmc_model = kmc.fit(traindata)
    return kmc_model


def KmeansClusterTest(testframe, model):
    testdata = np.array(testframe[CLUSTER_COL])
    kmc_result = model.predict(testdata)
    kmc_result = np.array(kmc_result).reshape(-1, 1)
    kmc_result = ['Class_' + str(int(kmc_result[i]) + 1)
                  for i in range(len(kmc_result))]
    testframe['cluster_result'] = kmc_result
    testframe = testframe.sort_values(by='cluster_result')
    return testframe


if __name__ == '__main__':
    fl = ['whole.csv',
          'stomach_si.csv',
          'prost_gland.csv',
          'skin.csv',
          'liver.csv',
          'bladder.csv',
          'colon.csv',
          'uters_uteri.csv',
          'breast.csv',
          'lung_bron.csv',
          'kidney.csv',
          'thyroid.csv',
          'ovary.csv']

    traindf = pd.read_csv(PATH + 'whole.csv')
    load_data(traindf)
    fig = plt.figure(figsize=(6, 5))
    for f in fl:
        testdf = pd.read_csv(PATH + f)
        load_data(testdf)
        testdf = testdf[CLUSTER_COL]
        for k in range(2, 10):
            h_model = KmeansClusterTrain(traindf, n_k=k)
            MEAND.append(sum(np.min(cdist(testdf,
                                          h_model.cluster_centers_,
                                          metric='euclidean'),
                                    axis=1))/testdf.shape[0])
        plt.plot(range(2, 10), MEAND)#, label=f.split('.')[0])
        plt.xlim(2, 9)
        plt.ylim(0, 0.5)
        MEAND = []
        print('==================')
    plt.tight_layout()
    plt.savefig(OUTPATH + 'bowl.png')
