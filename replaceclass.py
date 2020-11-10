import numpy as np
import pandas as pd
import os
import cfg
import argparse
CFG = cfg.CFG()

parser = argparse.ArgumentParser()
parser.add_argument('-f', type=str, help='the folder of files you want test.')
args = parser.parse_args()

OLDER_PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
    args.f + '/analysis_data_' + args.f + '/'


def ReplaceClass(datacsvname):
    dataframe = pd.read_csv(OLDER_PATH + datacsvname, dtype={'class': str})
    """
    dataframe['ICD_O3_pathology'] = dataframe['ICD_O3_pathology'].apply(
        lambda x: x[:3]).tolist()
        """
    classlist = list(set(dataframe['class']))
    for c in classlist:
        dataframe['class'][dataframe['class'] == c] = CFG.classreplace[c]
    print(list(set(dataframe['ICD_O3_pathology'])))
    dataframe.to_csv(OLDER_PATH + datacsvname, index=False)


if __name__ == '__main__':
    filelist = os.listdir(OLDER_PATH)
    for f in filelist:
        if 'csv' in f:
            print(f)
            ReplaceClass(f)
