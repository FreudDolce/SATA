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
from AMBscattersurvivalforgene import SurvivalAnalysis, SurvivalHeatMap

MEAND = []


CMAP = 'Blues'
N_CLUSTER = 7
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_patient_20200617/'
OUTPATH = '/Users/freud/Documents/MANU/GENE_ANALYSIS_1/manuscript/figures/Figure/f5/'
CLUSTER_COL = []
CLASS_LIST = []
for i in range(7):
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


def AnnovaTestAndDraw(annovaframe, test_col, title, class_col='cluster_result'):
    annovaframe.drop(annovaframe[annovaframe[test_col].isin(['data_loss', 'data_empty'])].index,
                     inplace=True)
    annovaframe[test_col] = pd.to_numeric(
        annovaframe[test_col], errors='coerce')
    class_list = list(set(annovaframe[class_col]))
    class_list = sorted(class_list)
    xLabel = []
    meanlist = []
    stdlist = []
    args = []
    for clas in class_list:
        xLabel.append(str(clas))
        classframe = annovaframe[annovaframe[class_col] == clas]
        args.append(classframe[test_col])
        classmean = classframe[test_col].mean()
        meanlist.append(classmean)
        classstd = classframe[test_col].std()
        stdlist.append(classstd)
    try:
        annova_f, annova_p = stats.f_oneway(*args)
    except ValueError:
        annova_f = 0
        annova_p = 'nan'
    fig = plt.figure(figsize=(len(class_list)/2 + 1, 4))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.bar(x=np.arange(len(class_list)),
           height=meanlist,
           width=0.8,
           yerr=stdlist,
           capsize=3)
    ax.tick_params(labelsize=30)
    print(xLabel)
    ax.set_xticks(range(len(xLabel)))
    ax.set_xticklabels(xLabel, rotation=90)
    ax.set_title(title + str(annova_p), fontsize=30)
    return (annova_f, annova_p, plt)


def DrawHeatMap(inframe):
    fig = plt.figure(dpi=300, figsize=(5, 4))
    Grid = plt.GridSpec(1, 30, wspace=0.1, hspace=0)
    ax_l = fig.add_subplot(Grid[0, 0])
    ax = fig.add_subplot(Grid[0, 1:])
    for c in CLASS_LIST:
        try:
            labelframe = pd.concat([labelframe, pd.DataFrame(np.ones((len(
                inframe[inframe['cluster_result'] == c]), 1)) * int(c.split('_')[1]))])
        except NameError:
            labelframe = pd.DataFrame(np.ones((len(
                inframe[inframe['cluster_result'] == c]), 1)) * int(c.split('_')[1]))
    ax_l.imshow(labelframe, cmap='Paired')
    ax_l.set_aspect('auto')
    ax_l.set_axis_off()
    plt.tight_layout()
    drawframe = inframe[CLUSTER_COL]
    im = ax.imshow(drawframe, cmap=CMAP, vmin=0, vmax=1, aspect='auto')
    ax.set_axis_off()
    # ax.set_xticks(range(0, 10))
    # ax.set_xticklabels(list(drawframe.columns), rotation=90)
    fig.colorbar(im)
    # plt.tight_layout
    # sns.heatmap(drawdf[CLUSTER_COL], vmin=0, vmax=1, cmap='OrRd')
    plt.tight_layout()
    return plt


if __name__ == '__main__':
    fl = ['stomach_si.csv',
          'whole.csv',
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
    for f in fl:
        if '.csv' in f:
            savepath = OUTPATH
            clusterdf = pd.read_csv(PATH + f)
            print(clusterdf)
            plt, survival = SurvivalAnalysis(
                clusterdf, classcol='cluster_result', sur_n=N_CLUSTER)
            plt.legend()
            plt.ylim(0, 1)
            plt.tick_params(labelsize=15)
            plt.xlim(0, )
            plt.xticks(rotation=90)
            plt.legend().set_visible(False)
            plt.tight_layout()
            plt.savefig(savepath + f.split('.')[0] + '_survival.tif')
            plt.close()
            plt = SurvivalHeatMap(survival, n=7, sig=0.99999)
            plt.tight_layout()
            plt.savefig(savepath + f.split('.')[0] + '_scatter.tif')
