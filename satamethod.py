# ===============================
# sata methods for patients
# Written by Ji Hongchen
# 20200909
# ===============================

import pandas as pd
import numpy as np
import os
import clincalcfg
from scipy import stats


PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp20200617/analysis_patient_20200617/'
CFG = clincalcfg.CFG()
CLASS_LIST = []
DROP_LIST = ['data_loss', 'data_empty']
STR_COL = ['patient_age', 'patient_weight']
for i in range(10):
    CLASS_LIST.append('Class_' + str(i + 1))


def load_data(dataframe):
    for c in STR_COL:
        dataframe[c] = dataframe[c].astype('str')
    for item in CFG.clicfeat_dict:
        for clas in CFG.clicfeat_dict[item]:
            dataframe[item][dataframe[item].isin(
                CFG.clicfeat_dict[item][clas])] = clas
    return dataframe


def Describe(dataframe, class_col):
    describe_item = ['count', 'max', 'min', 'q50%', 'mean', 'var', 'std']
    describe_list = []
    for cls in CLASS_LIST:
        for item in describe_item:
            describe_list.append(cls + '_' + item)
    describeframe = pd.DataFrame(columns=list(
        set(dataframe[class_col])), index=describe_list)
    for clicclass in describeframe.columns:
        _calcframe = dataframe[dataframe[class_col] == clicclass]
        for clas in CLASS_LIST:
            describeframe[clicclass].loc[clas +
                                         '_count'] = _calcframe['Rate'+clas].count()
            describeframe[clicclass].loc[clas +
                                         '_max'] = _calcframe['Rate'+clas].max()
            describeframe[clicclass].loc[clas +
                                         '_min'] = _calcframe['Rate'+clas].min()
            describeframe[clicclass].loc[clas +
                                         '_q50%'] = _calcframe['Rate'+clas].median()
            describeframe[clicclass].loc[clas +
                                         '_mean'] = _calcframe['Rate'+clas].mean()
            describeframe[clicclass].loc[clas +
                                         '_var'] = _calcframe['Rate'+clas].var()
            describeframe[clicclass].loc[clas +
                                         '_std'] = _calcframe['Rate'+clas].std()
    col_order = []
    for col in CFG.clicfeat_dict[class_col]:
        col_order.append(col)
        if col not in describeframe.columns:
            describeframe[col] = np.nan
    describeframe = describeframe[col_order]
    return describeframe


def AnnovaTest(dataframe, clac_item, class_list):
    dataframe = dataframe.drop(dataframe[dataframe[clac_item] == 'delete'].index)
    annovaframe = pd.DataFrame(columns=['F', 'P'], index=class_list)
    for clas in class_list:
        args = []
        clasframe = dataframe[[clac_item, 'Rate' + clas]]
        for clac_i in list(set(dataframe[clac_item])):
            args.append(clasframe[clasframe[clac_item]
                                  == clac_i]['Rate' + clas])
        f, p = stats.f_oneway(*args)
        annovaframe.loc[clas] = [f, p]
    return annovaframe

if __name__ == '__main__':
    df = pd.read_csv(
        '~/Documents/MANU/lstmsom_data/exp20200617/analysis_patient_20200617/whole.csv')
    ddf = load_data(df)
    dsf = AnnovaTest(ddf, clac_item='ajcc_stage', class_list=CLASS_LIST)
    print(dsf)
