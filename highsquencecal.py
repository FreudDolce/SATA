#=======================================
# Satistic genes has hige mut frequence
# Wrtten by Ji Hongchen
# 20200831
#=======================================

import pandas as pd
import numpy as np
import os

def CalFreqGene(dataframe):
    class_list = list(set(dataframe['class']))
    calcframe = pd.DataFrame()
    for i in range(10):
        classframe = dataframe[dataframe['class'] == 'Class_' + str(i+1)]
        freqframe = classframe['Hugo_Symbol'].value_counts()
        if len(freqframe) < 20:
            freqframe = freqframe.append(pd.Series([0] * (20 - len(freqframe))))
        print(freqframe)
        calcframe['Class_' + str(i+1)] = freqframe.index[0:20]
        calcframe['freq_' + str(i+1)] = list(freqframe[0:20])/freqframe.count()
    return calcframe

if __name__ == '__main__':
    file_list = os.listdir('/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/')
    for f in file_list:
        if '.csv' in f:
            df = pd.read_csv('~/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/' + f)
            cf = CalFreqGene(df)
            cf.to_csv('/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_data_20200617/' + f.split('.')[0] + '/high_freq.csv', index=False)

