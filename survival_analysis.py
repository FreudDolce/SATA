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
import os

GENEPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/other_files/'
DATAPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/'
OUTPUTPATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/genecalc/'
CLASS_LIST = []
LABEL_LIST = []
for i in range(10):
    CLASS_LIST.append('Class_' + str(i + 1))
    LABEL_LIST.append('Class ' + str(i + 1))


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
    plt.xticks(range(len(filelist)), filelist, rotation=90, fontsize=13)
    return plt


def SurvivalAnalysis(surframe, classcol, sur_n=10):
    kmf = KaplanMeierFitter()
    fig = plt.figure(figsize=(4.5, 4))
    survivalarray = np.ones((sur_n, sur_n))
    survivalframe = pd.DataFrame(survivalarray, columns=[
                                 CLASS_LIST], index=[CLASS_LIST])
    classlist = []
    for cla in CLASS_LIST:
        if cla in list(set(surframe[classcol])):
            classlist.append(cla)
    for i in range(len(classlist)):
        sataframe_1 = surframe[surframe[classcol] == classlist[i]]
        patient_num_1 = len(sataframe_1)
        #print('***********************8 n = ', patient_num_1)
        survivalmedian_i = sataframe_1['to_last_known_alive'].median()
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
            print('------------------------ m n = ', len(sataframe_2))
            calcframe_2 = sataframe_2.drop(sataframe_2[sataframe_2['Tumor_Sample_Barcode'].isin(
                list(set(sataframe_1['Tumor_Sample_Barcode'])))].index)
            print('++++++++++++++++++++++m n = ', len(calcframe_2))
            patient_num_2 = len(sataframe_2)
            survivalmedian_j = sataframe_2['to_last_known_alive'].median()
            _T_2 = calcframe_2['to_last_known_alive']
            _E_2 = calcframe_2['patient_status']
            lr_p = logrank_test(_T_1, _T_2, _E_1, _E_2, alpha=0.95)
            if lr_p.p_value != 1:
                if survivalmedian_j >= survivalmedian_i:
                    survivalframe.loc[classlist[j],
                                      classlist[i]] = -1 * lr_p.p_value
                    survivalframe.loc[classlist[i],
                                      classlist[j]] = lr_p.p_value
                else:
                    survivalframe.loc[classlist[j],
                                      classlist[i]] = lr_p.p_value
                    survivalframe.loc[classlist[i],
                                      classlist[j]] = -1 * lr_p.p_value
    return (plt, survivalframe)


def SurvivalScatter(survival):
    drawframe = survival.apply(abs)
    drawframe = drawframe.apply(np.log10)
    drawframe = drawframe.apply(abs)
    print(drawframe)
    fig = plt.figure(figsize=(4, 4))
    for i in range(len(drawframe)):
        for j in range(len(drawframe)):
            if 0 < survival.iloc[j, i] < 0.05:
                plt.scatter(
                    [i], [9 - j], s=[5 * drawframe.iloc[j, i] ** 2], c='crimson')
            elif 0 > survival.iloc[j, i] > -0.05:
                plt.scatter(
                    [i], [9 - j], s=[5 * drawframe.iloc[j, i] ** 2], c='steelblue')
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
        print ('>>>>>>>>', gene)
        plt = DrawPileMap(fl, gene)
        plt.tight_layout()
        plt.savefig(OUTPUTPATH + gene + '/DifferentInCancers.tif')
        print (gene + ', saved.')
    """
    geneframe = pd.read_csv(
        '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/other_files/whole.csv', index_col=0)
    calcgene = list(geneframe.index)[:20]
    for f in fl:
        dataframe = pd.read_csv(
            '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/' + f)
        load_data(dataframe)
        print(f)
        print(calcgene)
        print('=======================')
        for gene in calcgene:
            if os.path.exists(OUTPUTPATH + gene + '/' + f.split('.')[0]) == False:
                os.makedirs(OUTPUTPATH + gene + '/' + f.split('.')[0])
            df = dataframe[dataframe['Hugo_Symbol'] == gene]
            # print(len(df))
            #pl = df['Tumor_Sample_Barcode'].value_counts()
            # for p in pl.index:
            #    if pl[p] != 1:
            #        df.drop(df[df['Tumor_Sample_Barcode'] == p].index, inplace=True)
            # print(len(df))
            plt, survival = SurvivalAnalysis(df, classcol='class')
            plt.legend(bbox_to_anchor=[1.02, 1])
            plt.ylim(0, 1)
            plt.tick_params(labelsize=15)
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(OUTPUTPATH + gene + '/' +
                        f.split('.')[0] + '/lifeline.png')
            plt = SurvivalScatter(survival)
            plt.xlim(-0.5, 9.5)
            plt.ylim(-0.5, 9.5)
            plt.xticks(np.linspace(-0.5, 9.5, 11), [])
            plt.yticks(np.linspace(-0.5, 9.5, 11), [])
            plt.grid(linestyle='--', color='grey', alpha=0.4)
            # ax.spines['top'].set_visible(False)
            # ax.spines['right'].set_visible(False)
            # ax.spines['left'].set_visible(False)
            # ax.spines['bottom'].set_visible(False)
            #ax = plt.gca()
            # plt.axis('off')
            plt.tight_layout()
            plt.savefig(OUTPUTPATH + gene + '/' +
                        f.split('.')[0] + '/survivalana.png')

    """
