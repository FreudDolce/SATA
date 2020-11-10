# ================================
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

READAPTH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/'
SAVEPATH = '/Users/freud/Documents/MANU/GENE_ANALYSIS_1/manuscript/figures/'

AGE_ORDER = ['11 ~ 15', '16 ~ 20', '21 ~ 25', '26 ~ 30', '31 ~ 35',
             '36 ~ 40', '41 ~ 45', '46 ~ 50', '51 ~ 55', '56 ~ 60',
             '61 ~ 65', '66 ~ 70', '71 ~ 75', '76 ~ 80', '81 ~ 85', '86 ~ 90']

WEIGHT_ORDER = ['0 ~ 40', '41 ~ 50', '51 ~ 60', '61 ~ 70', '71 ~ 80',
                '81 ~ 90', '91 ~ 100', '101 ~ 110', '111 ~ 120', '>120', 'chi2', 'P', 'df']

AJCC_ORDER = ['I', 'II', 'III', 'IV']

CALC = 1000


def SortFrame(dataframe):
    class_list_costom = CFG.class_list_costom
    dataframe = dataframe.reindex(index=class_list_costom, fill_value=0)
    return dataframe


def DrawPileMap(dataframe, fign, rate='yes'):
    xLabel = list(dataframe.columns)
    #del xLabel[-3:]
    yLabel = list(dataframe.index)
    numframe = dataframe[xLabel]
    drawframe = dataframe[xLabel]
    fig = plt.figure(figsize=(len(xLabel)/2 + 1,
                              7))
    ax = fig.add_subplot(111)
    ax.set_ylim(0, 1)
    for i in drawframe.columns:
        if rate == 'yes':
            if drawframe[i].any() != 0:
                drawframe[i] = drawframe[i] / drawframe[i].sum()
    totalprecnt = [0] * len(xLabel)
    for i in yLabel:
        precnt = np.array(drawframe.loc[i])
        ax.bar(xLabel, precnt, width=0.95, bottom=totalprecnt, label=str(i))
        totalprecnt = totalprecnt + precnt
    for i in range(len(xLabel)):
        xLabel[i] += (' (n =' + str(int(dataframe[xLabel[i]].sum())) + ')')
        xLabel[i] = '0'
    plt.xticks(range(len(xLabel)), xLabel, rotation='vertical', size=35)
    plt.yticks([0, 0.25, 0.5, 0.75, 1], [
               '0%', '25%', '50%', '75%', '100%'], size=35)
    plt.tight_layout()
    savepath = SAVEPATH + fi + '/'
    if os.path.exists(savepath) == False:
        os.mkdir(savepath)
    plt.savefig(savepath + fign + '.tif')


if __name__ == '__main__':
    folder = ['whole.csv']#os.listdir(READAPTH)
    drawfile = 'ajcc_stage_base_as_item'
    fig_name = 'AJCC_STAGE'
    for fi in folder:
        if '.csv' not in fi:
            print(fi)
            try:
                drawfd = pd.read_csv(
                    READAPTH + fi + '/' + drawfile + '.csv', index_col=0)
                if 'weight' in drawfile:
                    order = WEIGHT_ORDER
                    for n in order:
                        try:
                            drawfd[n]
                        except KeyError:
                            drawfd[n] = 0
                    drawfd = drawfd[order]
                elif '_age' in drawfile:
                    order = AGE_ORDER
                    for n in order:
                        try:
                            drawfd[n]
                        except KeyError:
                            drawfd[n] = 0
                    drawfd = drawfd[order]
                elif 'ajcc' in drawfile:
                    order = AJCC_ORDER
                    for n in order:
                        try:
                            drawfd[n]
                        except KeyError:
                            drawfd[n] = 0
                    drawfd = drawfd[order]
                drawfd = SortFrame(drawfd)
                DrawPileMap(drawfd, fign=fig_name)

            except FileNotFoundError:
                print('FileNotFoundError')
                continue
