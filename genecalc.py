# ================================
# Sata mb of patients
# 20200807
# Written by Ji Hongchen
# ================================

import pandas as pd
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy import stats
import cfg
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-g', help="gene needed to be calc")
args = parser.parse_args()


CMAP = plt.get_cmap('Reds')
CFG = cfg.CFG()
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/'
CALC_PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/patients_mut/'
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
    dataframe.drop(dataframe[dataframe['patient_status'].isin(
        ['data_loss', 'data_empty'])].index, inplace=True)


def calc_mb(dataframe, f_name):
    patient_list = list(set(dataframe['Tumor_Sample_Barcode']))
    class_list = CLASS_LIST
    CalcFrame = pd.DataFrame(columns=class_list, index=patient_list)
    for patient in patient_list:
        patframe = dataframe[dataframe['Tumor_Sample_Barcode'] == patient]
        patmbcal = patframe['class'].value_counts()
        for clas in class_list:
            try:
                CalcFrame[clas].loc[patient] = patmbcal[clas] / \
                    (patframe['class'].count())
            except KeyError:
                CalcFrame[clas].loc[patient] = 0
    CalcFrame.to_csv(CALC_PATH + f_name + '.csv', index=False)
    return CalcFrame


def CalcAsGene(dataframe, gene, calcitem, classcol):
    calcframe = dataframe[dataframe['Hugo_Symbol'] == gene]
    if calcitem in ['patient_age', 'patient_weight']:
        fig = plt.figure(figsize=(3, 4), dpi=300)
        ax = fig.add_subplot(111)
        calcframe.drop(calcframe[calcframe[calcitem].isin(
            ['data_loss', 'data_empty'])].index, inplace=True)
        calcframe[calcitem] = pd.to_numeric(
            calcframe[calcitem], errors='coerce')
        classlist = list(set(calcframe[classcol]))
        mean = []
        std = []
        args = []
        count = []
        for cl in range(0, 10):
            clas = 'Class_' + str(cl + 1)
            if clas in classlist:
                count.append(
                    ' (n = ' + str(calcframe[calcframe[classcol] == clas][calcitem].count()) + ')')
                args.append(calcframe[calcframe[classcol] == clas][calcitem])
                mean.append(calcframe[calcframe[classcol]
                                      == clas][calcitem].mean())
                std.append(calcframe[calcframe[classcol]
                                     == clas][calcitem].std())
            else:

                mean.append(0)
                std.append(0)
                count.append(' (n = 0)')
        std = np.nan_to_num(std)
        try:
            annova_f, annova_p = stats.f_oneway(*args)
        except ValueError:
            annova_f, annova_p = (9999, 9999)
        ax.errorbar(x=np.arange(10),
                    y=mean,
                    yerr=std,
                    ecolor=('black'),
                    capsize=2,
                    fmt='co',
                    mfc='black',
                    mec='black')
        ax.set_title('P = ' + str(round(annova_p, 4)),
                     loc='right', fontsize=20)
        ax.set_xticks(range(10))
        LABEL = [LABEL_LIST[i] + count[i] for i in range(len(LABEL_LIST))]
        ax.set_xticklabels(labels=LABEL, rotation=90)
        ax.set_ylim(0,)
        plt.tick_params(labelsize=15)
        # ax.text(len(classlist), -1, 'p = ' + str(round(annova_p, 2)),
        #        ha='right')
    elif calcitem in ['patient_gender', 'ajcc_stage', 't_stage', 'n_stage', 'm_stage']:
        fig = plt.figure(figsize=(4, 4), dpi=300)
        ax = fig.add_subplot(111)
        calcframe.drop(calcframe[calcframe[calcitem].isin(
            CFG.clicfeat_dict[calcitem]['delete'])].index, inplace=True)
        for i in CFG.clicfeat_dict[calcitem]:
            calcframe[calcitem][calcframe[calcitem].isin(
                CFG.clicfeat_dict[calcitem][i])] = i
        crosstable = pd.crosstab(calcframe[calcitem],
                                 calcframe[classcol])
        cttable = np.array(crosstable)
        try:
            ctp = stats.chi2_contingency(cttable)[1]
        except ValueError:
            ctp = 9999
        bottom_list = [0] * len(CLASS_LIST)
        for clas in CLASS_LIST:
            if clas not in crosstable.columns:
                crosstable[clas] = 0
        crosstable = crosstable[CLASS_LIST]
        count = []
        for c in CLASS_LIST:
            count.append(' (n = ' + str(crosstable[c].sum()) + ')')
            crosstable[c] = crosstable[c]/crosstable[c].sum()
        for j in range(len(crosstable)):
            ax.bar(CLASS_LIST,
                   crosstable.iloc[j],
                   0.8,
                   bottom=bottom_list,
                   label=crosstable.index[j])
            bottom_list += crosstable.iloc[j]
        ax.set_title('P = ' + '%.03f' % ctp, loc='right', fontsize=20)
        ax.set_xticks(range(10))
        LABEL = [LABEL_LIST[i] + count[i] for i in range(len(LABEL_LIST))]
        ax.set_xticklabels(labels=LABEL, rotation=90)
        ax.legend(bbox_to_anchor=(1.02, 1))
        ax.set_ylim(0, 1)
        plt.tick_params(labelsize=15)
    plt.tight_layout()
    return plt




if __name__ == '__main__':
       for f in fl:
        df = pd.read_csv(PATH + f)
        for item in ['patient_age', 'patient_weight', 'ajcc_stage',
                     't_stage', 'n_stage', 'm_stage', 'patient_gender']:
            print(f)
            print(item)
            print(args.g)
            saveplt = CalcAsGene(
                df, gene=args.g, calcitem=item, classcol='class')
            savepath = CALC_PATH + args.g + '/' + f.split('.')[0] + '/'
            if os.path.exists(savepath) == False:
                os.makedirs(savepath)
            plt.savefig(savepath + item + '.png')
