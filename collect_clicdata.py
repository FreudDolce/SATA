# ================================================
# Used for collect data from different files
# Written by Ji Hongchen
# 20200329
# ================================================

import pandas as pd
import numpy as np
import os
import argparse

COLLECT_FILE = 'vote_logistic_regression.csv'
COLLECT_COLUMNS = ['lr_p']

parser = argparse.ArgumentParser()
parser.add_argument('-f', help="the folder input and ouput.")
args = parser.parse_args()

if args.f:
    PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
        args.f + '/analysis_data_' + args.f + '/'

if __name__ == '__main__':
    try:
        os.makedirs(PATH + 'collection')
    except FileExistsError:
        pass
    f_list = os.listdir(PATH)
    files_list = []
    for f in f_list:
        if '.csv' not in f:
            files_list.append(f)
    files_list.remove('collection')
    for item in COLLECT_COLUMNS:
        info_list = pd.read_csv(PATH + 'whole/' + COLLECT_FILE, index_col=0)
        resultframe = pd.DataFrame(columns=info_list.index, index=files_list)
        for files in files_list:
            info = pd.read_csv(PATH + files + '/' + COLLECT_FILE, index_col=0)
            for col in resultframe.columns:
                try:
                    resultframe[col][files] = info[item][col]
                except KeyError:
                    resultframe[col][files] = 'null'
        print(resultframe)
        resultframe.to_csv(PATH + 'collection/' +
                           COLLECT_FILE.split('.')[0] + '_' + item + '.csv')
