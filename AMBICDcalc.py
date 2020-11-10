#================================
# Used for draw figures
# Written by Ji Hongchen
# 20200727
# ================================
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import cfg
CFG = cfg.CFG()
CMAP = plt.get_cmap('Reds')
READAPTH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/'
SAVEPATH = '/Users/freud/Documents/MANU/GENE_ANALYSIS_1/manuscript/figures/Figure/f2/'
TYPICAL_GENE = 'TP53'
CURVE = False

AGE_ORDER = ['11 ~ 15', '16 ~ 20', '21 ~ 25', '26 ~ 30', '31 ~ 35',
             '36 ~ 40', '41 ~ 45', '46 ~ 50', '51 ~ 55', '56 ~ 60',
             '61 ~ 65', '66 ~ 70', '71 ~ 75', '76 ~ 80', '81 ~ 85', '86 ~ 90']

WEIGHT_ORDER = ['0 ~ 40', '41 ~ 50', '51 ~ 60', '61 ~ 70', '71 ~ 80',
                '81 ~ 90', '91 ~ 100', '101 ~ 110', '111 ~ 120', '>120', 'chi2', 'P', 'df']

AJCC_ORDER = ['I', 'II', 'III', 'IV']

CALC = 1000


def DrawPileMap(dataframe, figname, rate='yes'):
    xLabel = list(dataframe.columns)
    print(xLabel)
    dataframe.sort_index(inplace=True)
    yLabel = list(dataframe.index)
    drawframe = dataframe[xLabel]
    fig = plt.figure(figsize=(3.5, 6.5), dpi=300)
    ax = fig.add_subplot(111)
    ax.set_ylim(0, 1)
    for i in drawframe.columns:
        if rate == 'yes':
            if drawframe[i].any() != 0:
                drawframe[i] = drawframe[i] / drawframe[i].sum()
    totalprecnt = [0] * len(xLabel)
    print(drawframe)
    for i in yLabel:
        precnt = np.array(drawframe.loc[i])
        ax.bar(xLabel, precnt, width=0.95, bottom=totalprecnt, label=str(i))
        totalprecnt = totalprecnt + precnt
    for i in range(len(xLabel)):
        xLabel[i] += (' (n =' + str(int(dataframe[xLabel[i]].sum())) + ')')
    plt.xticks(range(len(xLabel)), xLabel, rotation='vertical', size=10)
    plt.yticks([0, 0.25, 0.5, 0.75, 1], [
               '0%', '25%', '50%', '75%', '100%'], size=10)
    plt.tight_layout()
    plt.savefig(SAVEPATH + figname + '.tif', dpi=300,
                format='tif', bbox_inches='tight')

if __name__ == '__main__':
    calcitem = 'ICD_O3_pathology'
    fl = ['whole.csv']#os.listdir(READAPTH)
    least_sata = 1000
    for f in fl:
        if '.csv' in f:
            filename = f.split('.')[0]
            df = pd.read_csv(READAPTH + f)
            df['class'].fillna('data_loss', inplace=True)
            df[calcitem].fillna('data_loss', inplace=True)
            df['class'][df['class'] == 'Class_10'] = 'Class_X'
            df.drop(df[df[calcitem].isin(
                ['data_loss', 'data_empty'])].index, inplace=True)
            for item in CFG.clicfeat_dict[calcitem]:
                df[calcitem][df[calcitem].isin(CFG.clicfeat_dict[calcitem][item])] = item
            df.drop(
                df[df['class'].isin(['data_loss', 'data_empty'])].index, inplace=True)
            crossframe = pd.crosstab(df['class'], df[calcitem])
            collist = list(crossframe.columns)
            collist.sort(reverse=True)
            crossframe = crossframe[collist]
            for col in crossframe.columns:
                if crossframe[col].sum() < least_sata:
                    crossframe.drop([col], axis=1, inplace=True)
            print(crossframe)
            DrawPileMap(crossframe, figname=filename+'_path_MBpile')
