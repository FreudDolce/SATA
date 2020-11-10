# ===============================
# Used for patient cluster
# 20200920
# Ji Hongchen
# ===============================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

GENEPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/other_files/'
DATAPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/'
OUTPUTPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/genecalc/'
CLASS_LIST = []
LABEL_LIST = []
for i in range(10):
    CLASS_LIST.append('Class_' + str(i + 1))
    LABEL_LIST.append('Class ' + str(i + 1))


def DrawPileMap(filelist, gene):
    indexlist = ['Class_' + str(i + 1) for i in range(10)]
    countframe = pd.DataFrame(columns=filelist, index=indexlist)
    for f in filelist:
        dataframe = pd.read_csv(DATAPATH + f)
        geneframe = dataframe[dataframe['Hugo_Symbol'] == gene]
        total_num = len(geneframe)
        MBseries = geneframe['class'].value_counts()
        for index in indexlist:
            try:
                countframe[f].loc[index] = MBseries[index] / total_num
            except KeyError:
                countframe[f].loc[index] = 0
    fig = plt.figure(figsize=(4, 5), dpi=300)
    bot = [0] * len(filelist)
    for index in indexlist:
        plt.bar(range(len(filelist)),
                countframe.loc[index],
                bottom=bot)
        bot += np.array(countframe.loc[index])
    plt.ylim(0, 1)
    plt.yticks([0, 0.25, 0.5, 0.75, 1], [0, 25, 50, 75, 100])
    plt.xticks(range(len(filelist)), filelist, rotation=90, fontsize=13)
    return plt


if __name__ == '__main__':
    genelist = ['TTN', 'MUC16', 'TP53',
                'SYNE1', 'CSMD3', 'LRP1B',
                'RYR2', 'DNAH5', 'OBSCN',
                'PCLO', 'USH2A', 'FAT4',
                'DST', 'ZFHX4', 'XIRP2']
    fl = ['whole.csv',
          'lung_bron.csv',
          'breast.csv',
          'prost_gland.csv',
          'colon.csv',
          'stomach_si.csv',
          'liver.csv',
          'uters_uteri.csv',
          'bladder.csv',
          'ovary.csv',
          'skin.csv',
          'kidney.csv',
          'thyroid.csv']
    for gene in genelist:
        print('>>>>>>>>', gene)
        plt = DrawPileMap(fl, gene)
        plt.tight_layout()
        plt.savefig(OUTPUTPATH + gene + '/DifferentInCancers.tif')
        print(gene + ', saved.')
