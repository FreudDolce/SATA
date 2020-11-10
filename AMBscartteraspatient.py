# =========================================
# Used for visualise the mutation sequence
# Written by Ji Hongchen
# 20200419
# =========================================

from cfg import CFG
import argparse
import os
from pyfaidx import Fasta
import SATA_PRETREAT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib


cfg = CFG()
parser = argparse.ArgumentParser()
parser.add_argument('-f', help="the folder input and ouput.")
args = parser.parse_args()

GENEDICT = {}
VISUALDICT = {}
EXT = 10  # cfg.ext
CAMP = matplotlib.cm.hot
WHOLE_SEQUENCE = Fasta(
    '/Users/freud/Documents/MANU/lstmsom_data/clinical_data/hg38.fa')
CALC_NUM = 3
CAL_NUM_CSV = 20
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/'
FIG_PATH = '/Users/freud/Documents/MANU/GENE_ANALYSIS_1/manuscript/figures/Figure/f2/'


def CountMBforPatients(dataframe, calcitem='ICD_O3_pathology'):
    class_list = ['Class_1', 'Class_2', 'Class_3', 'Class_4', 'Class_5',
                  'Class_6', 'Class_7', 'Class_8', 'Class_9', 'Class_10']
    claclist = list(set(dataframe[calcitem]))
    countframe = pd.DataFrame(columns=class_list)
    for i in claclist:
        if dataframe[dataframe[calcitem] == i][calcitem].count() > 1000:
            countframe.loc[i] = 0
    for clas in class_list:
        classframe = dataframe[dataframe['class'] == clas]
        for c in claclist:
            if c in countframe.index:
                try:
                    countframe[clas].loc[c] = \
                        classframe[classframe[calcitem] ==
                                   c]['Tumor_Sample_Barcode'].value_counts().mean()
                except KeyError:
                    countframe[clas].loc[c] = 0
    countframe.sort_index(inplace=True)
    return (countframe, countframe.index)


def RateMBforPatients(dataframe, ind, calcitem='ICD_O3_pathology'):
    class_list = ['Class_1', 'Class_2', 'Class_3', 'Class_4', 'Class_5',
                  'Class_6', 'Class_7', 'Class_8', 'Class_9', 'Class_10']
    rate_list = ['Rate' + class_list[i] for i in range(len(class_list))]
    rateframe = pd.DataFrame(columns=rate_list, index=ind)
    for k in cfg.clicfeat_dict[calcitem]:
        dataframe[calcitem][dataframe[calcitem].isin(
            cfg.clicfeat_dict[calcitem][k])] = k
    for cr in rateframe.index:
        crframe = dataframe[dataframe[calcitem] == cr]
        for i in rate_list:
            rateframe[i].loc[cr] = crframe[i].mean()
    return rateframe


def MBScatter(countframe, rateframe):
    print(countframe)
    print(rateframe)
    countframe = countframe.astype('int')
    countframe = countframe.apply(np.log2)
    fig = plt.figure(figsize=(10, 5), dpi=300)
    for i in range(len(countframe.columns)):
        for j in range(len(countframe.index)):
            c_1 = 0.7 + 0.3 * np.log10(rateframe.iloc[j, i])
            if c_1 < 0:
                c_1 = 0
            c_3 = -0.35 * np.log10(rateframe.iloc[j, i])
            if c_3 > 1:
                c_3 = 1
            plt.scatter([i], [j],
                        s=11.5 * (countframe.iloc[j, i] + 1),
                        c=(c_1, 0, c_3))
    yLabel = []
    for ind in countframe.index:
        yLabel.append((50 - len(ind)) * 'Z' + ind)
    plt.yticks(range(len(countframe.index)), yLabel, fontproperties='Courier New', size=13)
    return plt


if __name__ == '__main__':
    fl = os.listdir(PATH + 'analysis_data_20200617/')
    """
    countframe = np.array([1, 5, 10, 50, 100, 500])
    countframe = pd.DataFrame(countframe, columns=['a'], index = ['a', 'b', 'c', 'd', 'e', 'f'])
    print(countframe)
    rateframe = np.array([0.005, 0.01, 0.05, 0.1, 0.5])
    rateframe = pd.DataFrame(rateframe, columns=['a'], index = ['a', 'b', 'c', 'd','e'])
    print(rateframe)
    plt = MBScatter(countframe, rateframe)
    plt.tight_layout()
    plt.savefig(FIG_PATH + 'color_label.tif')
    """
    for f in fl:
        if '.csv' in f:
            print(f.split('.')[0])
            countf = PATH + 'analysis_data_20200617/' + f
            ratef = PATH + 'analysis_patient_20200617/' + f
            wholeframe = pd.read_csv(countf)
            ratedf = pd.read_csv(ratef)
            for p in cfg.clicfeat_dict['ICD_O3_pathology']:
                wholeframe['ICD_O3_pathology'].loc[wholeframe['ICD_O3_pathology'].isin(cfg.clicfeat_dict['ICD_O3_pathology'][p])] = p
            countframe, indexes = CountMBforPatients(
                wholeframe, calcitem='ICD_O3_pathology')
            rateframe = RateMBforPatients(
                ratedf, ind=indexes, calcitem='ICD_O3_pathology')
            plt = MBScatter(countframe, rateframe)
            plt.tight_layout()
            plt.savefig(FIG_PATH + f.split('.')[0] + '_p_scatter.tif')
