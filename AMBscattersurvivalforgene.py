# ===============================
# Used for patient cluster
# 20200920
# Ji Hongchen
# ===============================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import median_survival_times
import os

GENEPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/other_files/'
DATAPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/'
OUTPUTPATH = '/Users/freud/Documents/MANU/GENE_ANALYSIS_1/manuscript/figures/Figure/f3/'
CLASS_LIST = []
LABEL_LIST = []
for i in range(7):
    CLASS_LIST.append('Class_' + str(i + 1))
    LABEL_LIST.append('MB ' + str(i + 1))

LABELFRAME = pd.DataFrame(np.zeros((10, 10)))
for i in range(len(LABELFRAME)):
    LABELFRAME[i] = 5 - i


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


def SurvivalAnalysis(surframe, classcol, sur_n=10):
    kmf = KaplanMeierFitter()
    fig = plt.figure(figsize=(4, 4), dpi=300)
    survivalarray = np.ones((sur_n, sur_n))
    survivalframe = pd.DataFrame(survivalarray, columns=[
                                 CLASS_LIST], index=[CLASS_LIST])
    classlist = []
    for cla in CLASS_LIST:
        if cla in list(set(surframe[classcol])):
            classlist.append(cla)
    classlist = np.array(classlist)
    for i in range(len(classlist)):
        sataframe_1 = surframe[surframe[classcol] == classlist[i]]
        patient_num_1 = len(sataframe_1)
        _T_1 = sataframe_1['to_last_known_alive']
        _E_1 = sataframe_1['patient_status']
        kmf.fit(_T_1, event_observed=_E_1,
                label=' '.join(classlist[i].split('_')) + '   (n = ' + str(patient_num_1) + ')')
        kmf.plot(ci_alpha=0)
        for j in range(i+1, len(classlist)):
            sataframe_2 = surframe[surframe[classcol] == classlist[j]]
            calcframe_1 = sataframe_1.drop(sataframe_1[sataframe_1['Tumor_Sample_Barcode'].isin(
                list(set(sataframe_2['Tumor_Sample_Barcode'])))].index)
            _T_1 = calcframe_1['to_last_known_alive']
            _E_1 = calcframe_1['patient_status']
            try:
                kmf.fit(_T_1, event_observed=_E_1)
                survivalmedian_i = kmf.median_survival_time_
            except ValueError:
                survivalmedian_i = 0
            calcframe_2 = sataframe_2.drop(sataframe_2[sataframe_2['Tumor_Sample_Barcode'].isin(
                list(set(sataframe_1['Tumor_Sample_Barcode'])))].index)
            patient_num_2 = len(sataframe_2)
            _T_2 = calcframe_2['to_last_known_alive']
            _E_2 = calcframe_2['patient_status']
            try:
                kmf.fit(_T_2, event_observed=_E_2)
                survivalmedian_j = kmf.median_survival_time_
            except ValueError:
                survivalmedian_j = 0
            try:
                lr_p = logrank_test(_T_1, _T_2, _E_1, _E_2)
            except UnboundLocalError:
                lr_p.p_value = 1
            if lr_p.p_value != 1:
                if survivalmedian_j >= survivalmedian_i:
                    survivalframe.iloc[j, i] = -1 * lr_p.p_value
                    survivalframe.iloc[i, j] = lr_p.p_value
                else:
                    survivalframe.iloc[j, i] = lr_p.p_value
                    survivalframe.iloc[i, j] = -1 * lr_p.p_value
    return (plt, survivalframe)


def SurvivalHeatMap(survival, n=10, sig=0.05):
    for i in survival.columns:
        for j in survival.index:
            if ((survival[i][j] >= sig) or (survival[i][j] <= -sig)):
                survival[i][j] = 0
            elif survival[i][j] > 0:
                survival[i][j] = np.log10(survival[i][j])
            elif survival[i][j] <= 0:
                survival[i][j] = -1 * np.log10(-1 * survival[i][j])
    survival.fillna(0)
    for i in survival.columns:
        for j in survival.index:
            if survival[i][j] > 4:
                survival[i][j] = 4
            elif survival[i][j] < -4:
                survival[i][j] = -4
    fig = plt.figure(figsize=(4, 4), dpi=300)
    plt.xlim(-0.5, n-0.5)
    plt.ylim(-0.5, n-0.5)
    plt.grid(color='grey', alpha=0.4)
    plt.xticks(np.linspace(-0.5, n-0.5, n+1), [])
    plt.yticks(np.linspace(-0.5, n-0.5, n+1), [])
    plt.imshow(survival, cmap='RdBu', vmin=-4, vmax=4)
    return plt


def SurvivalScatter(survival):
    drawframe = survival.apply(abs)
    drawframe = drawframe.apply(np.log10)
    drawframe = drawframe.apply(abs)
    fig = plt.figure(figsize=(4, 4))
    for i in range(len(drawframe)):
        for j in range(len(drawframe)):
            if 0 < survival.iloc[j, i]: # <= 0.05:
                plt.scatter(
                    [i], [9 - j], s=[5 * drawframe.iloc[j, i] ** 2], c='orangered')
            elif 0 > survival.iloc[j, i]: # > -0.05:
                plt.scatter(
                    [i], [9 - j], s=[5 * drawframe.iloc[j, i] ** 2], c='steelblue')
    return plt


if __name__ == '__main__':
    calcgene = ['TTN', 'MUC16', 'TP53',
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
    for f in fl:
        print('------------------------>', f)
        dataframe = pd.read_csv(
            '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/' + f)
        load_data(dataframe)
        print('=======================')
        for gene in calcgene:
            print(gene)
            print('=====================================')
            if os.path.exists(OUTPUTPATH + gene + '/') == False:
                os.makedirs(OUTPUTPATH + gene + '/')
            df = dataframe[dataframe['Hugo_Symbol'] == gene]
            print(len(df))
            pl = df['Tumor_Sample_Barcode'].value_counts()
            # for p in pl.index:
            #    if pl[p] != 1:
            #        df.drop(df[df['Tumor_Sample_Barcode']
            #                   == p].index, inplace=True)
            print(len(df))
            _, survival = SurvivalAnalysis(df, classcol='class')
            #plt.legend(bbox_to_anchor=(1, 1.4))
            #plt.tight_layout()
            #plt.savefig(OUTPUTPATH + 'label.tif', bbox_inches='tight')
            plt = SurvivalHeatMap(survival)
            plt.tight_layout()
            plt.savefig(OUTPUTPATH + gene + '/' + f.split('.')
                        [0] + '_survivalheatmap.tif')

