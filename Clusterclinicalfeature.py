# ==================================
# Used for draw figures by seaborn
# Ji Hongchen
# 20200911
# ==================================

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import satamethod
import shutil
import cfg
from scipy import stats

CFG = cfg.CFG()
N_CLUSTER = 7
CLASS_LIST = []
for i in range(N_CLUSTER):
    CLASS_LIST.append('Class_' + str(i + 1))
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_patient_20200617/'
CLAC_DICT = {'patient_age': 'avr', 'patient_weight': 'avr', 'patient_gender': 'chi',
             'ajcc_stage': 'chi', 't_stage': 'chi', 'n_stage': 'chi', 'm_stage': 'chi'}
OUTPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/cluster/'


def data_loader(dataframe, col_item):
    dataframe.drop(dataframe[dataframe[col_item].isin(
        CFG.clicfeat_dict[col_item]['delete'])].index, inplace=True)
    if CLAC_DICT[col_item] == 'chi':
        for i in CFG.clicfeat_dict[col_item]:
            dataframe[col_item][dataframe[col_item].isin(
                CFG.clicfeat_dict[col_item][i])] = i
    return dataframe


def HeatMap(dataframe, col_item):
    crosstab = pd.crosstab(dataframe[col_item], dataframe['cluster_result'])
    classlist = crosstab.columns
    try:
        chi_2, p, free_n = stats.chi2_contingency(np.array(crosstab))[0:3]
    except ValueError:
        chi_2 = 0
        p = 'nan'
        free_n = 0
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    xLabel = crosstab.columns
    for i in xLabel:
        crosstab[i] = crosstab[i]/crosstab[i].sum()
    bot = [0.0] * len(crosstab.columns)
    for item in crosstab.index:
        ax.bar(xLabel, list(crosstab.loc[item]),
               bottom=bot)
        bot = [bot[j] + crosstab.loc[item][j]
               for j in range(len(crosstab.columns))]
    ax.set_xticks(range(len(classlist)))
    ax.set_xticklabels(classlist, rotation=90)
    ax.legend(bbox_to_anchor=[1.02, 1])
    ax.set_title('p = ' + str(p), fontsize=20)
    return p, plt


def AnnovaTestAndDraw(annovaframe, test_col, class_col='cluster_result'):
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
    fig = plt.figure(figsize=(4,5))
    ax = fig.add_subplot(111)
    ax.bar(x=np.arange(len(class_list)),
           height=meanlist,
           width=0.8,
           yerr=stdlist,
           capsize=3)
    ax.tick_params(labelsize=10)
    print(xLabel)
    ax.set_xticks(range(len(xLabel)))
    ax.set_xticklabels(xLabel, rotation=90)
    ax.set_title(str(annova_p), fontsize=30)
    return (annova_f, annova_p, plt)


if __name__ == '__main__':
    fl = ['whole.csv']
    for f in fl:
        savepath = OUTPATH
        testdf = pd.read_csv(PATH + f)
        for calcitem in CLAC_DICT:
            testdf = data_loader(testdf, calcitem)
            if CLAC_DICT[calcitem] == 'chi':
                p, plt = HeatMap(testdf, calcitem)
                plt.ylim(0, 1)
            elif CLAC_DICT[calcitem] == 'avr':
                a, p, plt = AnnovaTestAndDraw(testdf,
                                              test_col=calcitem)
            plt.savefig(savepath + f.split('.')
                        [0] + '_' + calcitem + '.tif', bbox_inches='tight')
            plt.close()
