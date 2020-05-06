# ====================================
# Used for satastic
# Written by Ji Hongchen
# 20200216
# ====================================

import pandas as pd
import numpy as np
import os
#from sata_main import PARA


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


def Classifier(dataframe, accordcol, thresh=[0]):
    """
    return a dataframe with a column named "class_result"
    Should enter accordcol (the column according to )
    and thres (the list of threshes, upper line should be contained)
    """
    threshold = sorted(thresh)[:: -1]
    dataframe['class_result'] = -2
    for i in range(len(threshold)):
        dataframe['class_result'][dataframe[accordcol]
                                  <= threshold[i]] = len(threshold) - i
    dataframe['class_result'][dataframe['class_result'] == -2] = 0
    print('===================================================================')
    print('Hint: add a new column in dataframe, named "class_result"')
    print('===================================================================')
    return dataframe


def ClassiVoter(dataframe, freq_dict={}, method='quadrant', bias=[-0.2, 0.2],
                vote_line=1.2, direct=np.array([[0, 1], [1, 0], [0, -1], [-1, 0]])):
    print('******-> This frame should contain "class_result" column <-********')
    dataframe['vote_class'] = None
    id_list = np.unique(np.array(dataframe['Tumor_Sample_Barcode']))
    calc_dict = {}
    for patient_id in id_list:
        vote_freq = dataframe['class_result'][dataframe['Tumor_Sample_Barcode'] ==
                                              patient_id].value_counts()
        vote_freq = vote_freq/vote_freq.sum()
        try:
            calcframe.loc[patient_id] = 0
        except UnboundLocalError:
            calcframe = pd.DataFrame(columns=list(freq_dict))
            calcframe.loc[patient_id] = 0
        for i in vote_freq.index:
            calc_dict[i] = vote_freq[i]/freq_dict[i]
            calcframe[i][patient_id] = calc_dict[i]
        if method == 'fold':
            max_class = max(calc_dict, key=calc_dict.get)
            if calc_dict[max_class] >= vote_line:
                vote_class = max_class
            else:
                vote_class = 0
        elif method == 'quadrant':
            judgemati = [0, 0]
            for i in range(len(direct)):
                try:
                    judgemati += direct[i]*calc_dict[i + 1]
                except ValueError:
                    print(direct)
                    print(calc_dict)
                except KeyError:
                    judgemati += [0, 0]
            if (judgemati[0] >= vote_line + bias[0]) and (judgemati[1] > vote_line + bias[1]):
                vote_class = '1'
            elif (judgemati[0] > vote_line + bias[0]) and (judgemati[1] <= -1 * vote_line + bias[1]):
                vote_class = '2'
            elif (judgemati[0] <= -1 * vote_line + bias[0]) and (judgemati[1] < -1 * vote_line + bias[1]):
                vote_class = '3'
            elif (judgemati[0] < -1 * vote_line + bias[0]) and (judgemati[1] >= vote_line + bias[1]):
                vote_class = '4'
            else:
                vote_class = '0'
        dataframe['vote_class'][dataframe['Tumor_Sample_Barcode'] ==
                                patient_id] = vote_class
        clinical_series = dataframe[dataframe['Tumor_Sample_Barcode'] == patient_id]
        try:
            vote_frame = vote_frame.append(clinical_series.iloc[[0]])
        except UnboundLocalError:
            vote_frame = clinical_series.iloc[[0]]
    calcframe.fillna(0)
    print('******-> Important!!! retun a new dataframe dif from ori <-********')
    print('===================================================================')
    print('Hint: add a new column in dataframe, named "vote_class"')
    print('===================================================================')
    return (vote_frame, calcframe)


def DropEmptyData(dataframe, column, empty_index=['data_loss', 'data_empty']):
    for word in empty_index:
        dataframe = dataframe.drop(dataframe[dataframe[column] == word].index)
    return dataframe


def NeuPreteate(dataframe):
    neumeric_list = ['pred_1', 'pred_2', 'pred_3', 'pred_4',
                     'pred_5', 'pred_6', 'pred_7', 'pred_8']
    for item in neumeric_list:
        dataframe[item] = pd.to_numeric(dataframe[item], errors='coerce')
    return dataframe


def ChiSquarePretreat(dataframe, column, class_column='class_result'):
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
    pass
    C_22_ori = pd.read_csv(
        '~/Documents/MANU/lstmsom_data/exp20200227/cluster_result_20200227/flag_C02.csv')
    C_22_cla = Classifier(C_22_ori)
    C_22_cla = ClassiVoter(C_22_cla)
    print(C_22_cla)
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
