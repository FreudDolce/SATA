# ==================================
# Used for draw figures
# Ji Hongchen
# 20200911
# ==================================

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import satamethod
import shutil
import clincalcfg

CFG = clincalcfg.CFG()
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_patient_20200617/'
OUTPATH = '/Users/freud/Documents/MANU/GENE_ANALYSIS_1/manuscript/figures/Figure/f4/'
CMAP = plt.get_cmap('RdBu')

HEATLIST = []
CLASSLIST = []
for i in range(10):
    CLASSLIST.append('Class_' + str(i + 1))
    HEATLIST.append('RateClass_' + str(i + 1))


def loadcalcudata(dataframe, extract=['mean', 'std']):
    if len(extract) != 0:
        for na in dataframe.index:
            if na.split('_')[-1] not in extract:
                dataframe.drop([na], inplace=True)
    del dataframe['delete']
    return dataframe


def DrawCruveMap(oriframe, sataitem, fit=False):
    dataframe = satamethod.Describe(oriframe, class_col=sata_item)
    dataframe = loadcalcudata(dataframe)
    # for cn in list(CFG.clicfeat_dict[sataitem]):
    #    if cn not in dataframe.columns:
    #        dataframe[cn] = 0
    fig = plt.figure(figsize=(len(CFG.clicfeat_dict[sataitem])/2 + 1,
                              40), dpi=200)
    drawlist = dataframe.index
    for clas in CLASSLIST:
        ax = fig.add_subplot(10, 1, int(clas.split('_')[1]))
        countseries = oriframe[sataitem].dropna().value_counts()
        xlabel = []
        for name in dataframe.columns:
            try:
                xlabel.append(
                    str(name) + ' (n = ' + str(countseries[name]) + ')')
            except KeyError:
                xlabel.append(str(name) + ' (n = 0)')
        ax.set_xticks(range(len(dataframe.columns)))

        dataframe = dataframe.fillna(-1)
        # if clas.split('_')[1] in ['1', '2', '6', '9', '10']:
        #    ax.set_ylim(0, 0.15)
        # elif clas.split('_')[1] in ['3', '5', '7']:
        #    ax.set_ylim(0, 0.3)
        # else:
        #    ax.set_ylim(0, 0.6)
        bar_x = np.arange(len(dataframe.columns))
        y = dataframe.loc[clas+'_mean']
        y_top = np.round(y.max(), decimals=1)
        ax.errorbar(x=bar_x, y=y,
                    yerr=dataframe.loc[clas+'_std'],
                    fmt='co',
                    marker='o',
                    mfc='w',
                    mec='k',
                    ms=15,
                    mew=4,
                    ecolor='k',
                    elinewidth=4)
        ax.tick_params(labelsize=30)
        try:
            ax.set_yticks(np.linspace(0, y_top + 0.1, int(10 * y_top + 2)))
        except ValueError:
            ax.set_yticks(np.linspace(0, 1, 11),)
        try:
            ax.set_yticklabels(np.round(np.linspace(
                0, y_top + 0.1, int(10 * y_top + 2)), 1), fontsize=40)
        except ValueError:
            ax.set_yticklabels(np.round(np.linspace(0, 1, 11), 1), fontsize=40)
        ax.set_ylim(bottom=0, top=y_top+0.1)
        ax.set_xlim(-0.5, len(dataframe.columns) - 0.5)
        if clas == 'Class_10':
            ax.set_xticklabels(xlabel, rotation=90, fontsize=40,)
        else:
            ax.set_xticklabels([])
        if fit == True:
            x_fit = list(range(len(dataframe.columns)))
            x_layout = np.arange(-0.5, len(dataframe.columns), 0.1)
            y_fit = list(y)
            f = np.polyfit(x_fit, y_fit, 1)
            p = np.poly1d(f)
            y_layout = p(x_layout)
            ax.plot(x_lapyout, y_layout, 'r')
            plt.tight_layout()
    return plt


if __name__ == '__main__':
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
    # for f in fl:
    #    if '.csv' not in f:
    #        shutil.rmtree(PATH + f)
    # fl = ['breast.csv'] #os.listdir(PATH)
    # , 'patient_age', 'patient_gender', 't_stage', 'patient_weight', 'ajcc_stage', 'n_stage', 'm_stage']:
    for f in fl:
        if '.csv' in f:
            df = pd.read_csv(PATH + f)
            for sata_item in ['t_stage', 'patient_age', 'patient_gender', 't_stage', 'patient_weight', 'ajcc_stage', 'n_stage', 'm_stage']:
                df = satamethod.load_data(df)
                print('----->', sata_item)
                try:
                    annovaframe = satamethod.AnnovaTest(df,
                                                        clac_item=sata_item,
                                                        class_list=CLASSLIST)
                except TypeError:
                    annovaframe = pd.DataFrame([[0, 0], [0, 0]])
                print(f)
                print('===============================================================')
                # PatientHeatMap(df,
                #               output_path=PATH +
                #               f.split('.')[0] + '/heat_' + sata_item + '.tif',
                #               chose_col=HEATLIST,
                #               calc_item=sata_item)
                if os.path.exists(OUTPATH + sata_item) == False:
                    os.mkdir(OUTPATH + sata_item)
                annovaframe.to_csv(OUTPATH + sata_item + '/' +f.split('.')
                                   [0] + '_annova.csv')
                plt = DrawCruveMap(oriframe=df,
                                   fit=False,
                                   sataitem=sata_item)
                plt.tight_layout()
                plt.savefig(OUTPATH + sata_item + '/' +f.split('.') [0] + '_tif')
