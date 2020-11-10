# =======================================
# Key gene satastic
# Written by Ji Hongchen
# 20200704
# =======================================

import pandas as pd
import numpy as np
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300

parser = argparse.ArgumentParser()
parser.add_argument('-f', help="the folder input and ouput.")
args = parser.parse_args()
CLASS = ['Class_' + str(i + 1) for i in range(10)]
print(CLASS)

PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
    str(args.f) + '/analysis_data_' + str(args.f) + '/'
NUM_CAL = 30
CMAP = 'OrRd'
NORM = matplotlib.colors.Normalize(vmin=0, vmax=1)


def CalGeneNum(mutfile):
    mutinfo = pd.read_csv(mutfile)
    mutinfo = mutinfo[['Hugo_Symbol', 'class']]
    claslist = list(set(np.array(mutinfo['class'])))
    gene_count = mutinfo['Hugo_Symbol'].value_counts()
    tpklgenelist = list(gene_count[0:NUM_CAL].index)
    tpklgeneframe = pd.DataFrame(columns=claslist, index=tpklgenelist)
    for n in claslist:
        choseeninfo = mutinfo[mutinfo['class'] == n]
        for i in tpklgenelist:
            try:
                tpklgeneframe[n][i] = int(
                    choseeninfo['Hugo_Symbol'].value_counts()[i])
            except KeyError:
                tpklgeneframe[n][i] = 0
    tpklgeneframe.loc['Col_Sum'] = tpklgeneframe.apply(lambda x: x.sum())
    for n in claslist:
        try:
            tpklgeneframe[n] = tpklgeneframe[n] / tpklgeneframe[n]['Col_Sum']
        except ZeroDivisionError:
            tpklgeneframe[n] = 0
    tpklgeneframe.drop(['Col_Sum'], inplace=True)
    return tpklgeneframe


def DrawHeatMap(dataframe, savepath):
    dataframe = dataframe[CLASS]
    xlabel = dataframe.columns
    ylabel = dataframe.index
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yticks(range(len(ylabel)))
    ax.set_yticklabels(ylabel)
    ax.set_xticks(range(len(xlabel)))
    ax.set_xticklabels(xlabel)
    im = ax.imshow(np.array(dataframe).astype('float32'), cmap=CMAP)
    plt.colorbar(im)
    plt.savefig(savepath, dpi=300)


if __name__ == '__main__':
    file_list = os.listdir(PATH)
    for fi in file_list:
        if '.csv' in fi:
            print('{:=^120}'.format(fi + '_typical_gene_sata'))
            tpklgeneframe = CalGeneNum(PATH + fi)
            tpklgeneframe.to_csv(PATH + fi.split('.')[0] + '/genefreq.csv')
            DrawHeatMap(dataframe=tpklgeneframe,
                        savepath=PATH + fi.split('.')[0] + '/genefreq.png')
