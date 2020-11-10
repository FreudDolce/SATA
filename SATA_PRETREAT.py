# ====================================
# Used for satastic
# Written by Ji Hongchen
# 20200216
# ====================================

import pandas as pd
import numpy as np
import os
import cfg
#from sata_main import PARA

CFG = cfg.CFG()


def GetSurvTime(dataframe):
    """
    return a data frame of the same shape with input data
    to_last_known_alive is the survival tiem
    last_tumor_statu is the statu of death (0:live, 1:dead)
    """
    dataframe['to_last_known_alive'][dataframe['to_last_followup'] !=
                                     'data_loss'] = dataframe['to_last_followup']
    dataframe['to_last_known_alive'][dataframe['to_death_time'] !=
                                     'data_loss'] = dataframe['to_death_time']
    dataframe['to_last_known_alive'][dataframe['to_death_time'] ==
                                     'data_empty'] = dataframe['to_last_followup']
    dataframe['to_last_known_alive'][dataframe['to_last_followup'] ==
                                     'data_empty'] = dataframe['to_death_time']
    dataframe['patient_status'][dataframe['to_death_time'] !=
                                'data_loss'] = 1
    dataframe['patient_status'][dataframe['to_death_time'] ==
                                'data_loss'] = 0
    print('===================================================================')
    print('death statu is in column "last_tumor_statu"')
    print('Use survival time in "to_last_known_alive"')
    print('===================================================================')
    return dataframe


def Classifier(dataframe):
    """
    return a dataframe with a column named "class_result"
    Should enter accordcol (the column according to )
    and thres (the list of threshes, upper line should be contained)
    """
    dataframe.rename(columns={'class': 'class_result'}, inplace=True)
    print('===================================================================')
    print('Hint: add a new column in dataframe, named "class_result"')
    print('===================================================================')
    return dataframe


def ClassiVoter(dataframe, class_cal=122, line=2):
    print('******-> Important!!! retun a new dataframe dif from ori <-********')
    dataframe = dataframe[CFG.clinical_item]
    patient_list = list(set(np.array(dataframe['Tumor_Sample_Barcode'])))
    for patient in patient_list:
        try:
            judge_count = dataframe[dataframe['Tumor_Sample_Barcode']
                                    == patient]['class_result'].value_counts()[class_cal]
        except KeyError:
            judge_count = 0
        if judge_count >= line:
            dataframe['class_result'][dataframe['Tumor_Sample_Barcode']
                                      == patient] = str(class_cal) + '>=' + str(line)
        else:
            dataframe['class_result'][dataframe['Tumor_Sample_Barcode']
                                      == patient] = 'negtive'
    cli_frame = dataframe.drop_duplicates(subset='Tumor_Sample_Barcode')
    cli_frame['to_last_known_alive'] = pd.to_numeric(
        cli_frame['to_last_known_alive'], errors='coerce')
    cli_frame.drop(
        cli_frame[cli_frame['to_last_known_alive'] < CFG.min_alive].index, inplace=True)
    print('===================================================================')
    print('Hint: add a new column in dataframe, named "vote_class"')
    print('===================================================================')
    return cli_frame


def DropEmptyData(dataframe, column, empty_index=['data_loss', 'data_empty']):
    for word in empty_index:
        dataframe = dataframe.drop(dataframe[dataframe[column] == word].index)
    return dataframe


def NeuPreteate(dataframe, neumeric_list = ['pred_1', 'pred_2', 'pred_3', 'pred_4',
                     'pred_5', 'pred_6', 'pred_7', 'pred_8', 'ICD_O3_pathology']):
    for item in neumeric_list:
        dataframe[item] = pd.to_numeric(dataframe[item], errors='coerce')
    return dataframe


def ChiSquarePretreat(dataframe, column, class_column):
    class_list = np.unique(np.array(dataframe[class_column]))
    for i in class_list:
        part_frame = dataframe[dataframe[class_column] == i]
        count_series = part_frame[column].value_counts()
        count_frame = pd.DataFrame(count_series.values.reshape(1, -1),
                                   columns=count_series.index)
        count_frame['class_result'] = i
        try:
            chisquareframe = chisquareframe.append(count_frame)
        except UnboundLocalError:
            chisquareframe = count_frame
    chisquareframe.index = np.array(chisquareframe['class_result'])
    chisquareframe = chisquareframe.T
    chisquareframe.drop(['class_result'], inplace=True)
    return chisquareframe


if __name__ == '__main__':
    C_22_ori = pd.read_csv(
        '/Users/freud/Documents/MANU/lstmsom_data/exp20200616/analysis_data_20200616/adrenal.csv')
    C_22_cla = Classifier(C_22_ori)
    C_22_vo = ClassiVoter(C_22_cla)
    print(C_22_vo)
    """
    files = os.listdir(
        '/Users/freud/Documents/MANU/lstmsom_data/cluster_result_20200215/')
    for f in files:
        print(f, ' loading...')
        df = pd.read_csv('/Users/freud/Documents/MANU/lstmsom_data/cluster_result_20200215/'+f,
                    converters={'to_death_time': str, 'to_last_followup': str,
                                'to_last_known_alive': str})
        df_f = CalSurvTime(df)
        df_f.to_csv('/Users/freud/Documents/MANU/lstmsom_data/sata_data_20200215/'+f,
                    index=False)
        print(f, ' loaded')
    """
