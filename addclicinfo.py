#===================================
# Used for giving cluster result
# Written by Ji Hongchen
# 20200213
# ===================================

import numpy as np
import pandas as pd
import os
import GETPERDATA
import json
import argparse
from datetime import datetime
import time
import cfg

CFG = cfg.CFG()

COMBINE_DICT = {'head_neck':    ['C00', 'C01', 'C02', 'C03', 'C04',
                                 'C05', 'C06', 'C07', 'C08', 'C09',
                                 'C10', 'C13', 'C14'],
                'esophagus':    ['C15'],
                'stomach_si':   ['C16', 'C17'],
                'colon':        ['C18'],
                'rectom_anus':  ['C19', 'C20', 'C21'],
                'liver':        ['C22'],
                'baliary':      ['C24'],
                'pancerase':    ['C25'],
                'nasal':        ['C30', 'C32'],
                'lung_bron':    ['C34'],
                'thymus':       ['C37'],
                'heart_pleura': ['C38'],
                'bone_joints':  ['C40', 'C41'],
                'skin':         ['C44'],
                'peritoneum':   ['C48'],
                'breast':       ['C50'],
                'girl_genital': ['C51', 'C52', 'C53'],
                'uters_uteri':  ['C54', 'C55'],
                'ovary':        ['C56'],
                'prost_gland':  ['C61'],
                'testis':       ['C62'],
                'kidney':       ['C64'],
                'bladder':      ['C67'],
                'eyes_adnexa':  ['C69'],
                'meni_brain':   ['C70', 'C71'],
                'thyroid':      ['C73'],
                'adrenal':      ['C74'],
                'lymph_node':   ['C77']}

parser = argparse.ArgumentParser()
parser.add_argument('-f', type=str, help='the folder of files you want test.')
args = parser.parse_args()

args.f = str(args.f)

INFO_FROM = '/home/ji/Documents/test_data/exp' + args.f + '/mut_info/'
OLDER_PATH = '/home/ji/Documents/test_data/exp' + \
    args.f + '/analysis_data_' + args.f + '/'
CLINICAL_INFO_PATH = '/home/ji/Documents/lstmsom_data/clinical_data/clinical_info.csv'
FILE_LIST = np.load(
   '/home/ji/Documents/lstmsom_data/clinical_data/mutfilelist.npy')
CLASS_COL = ['ICD_O3_site']

cmd_1 = 'cat ' + INFO_FROM + '*.csv > ' + INFO_FROM + 'whole.csv'
os.system(cmd_1)

def GetClinclInfo(path):
    clinic_info = pd.read_csv(path)
    return (clinic_info)


if __name__ == '__main__':
    clinic_info = GetClinclInfo(CLINICAL_INFO_PATH)
    clinic_col = clinic_info.columns
    whole_info = pd.read_csv(INFO_FROM + 'whole.csv')
    whole_info.drop(whole_info[whole_info['Tumor_Sample_Barcode'] == 'Tumor_Sample_Barcode'].index, inplace=True)
    for c_i in clinic_col:
        whole_info[c_i] = 0
    patient_list = list(set(np.array(whole_info['Tumor_Sample_Barcode'])))
    num = 0
    for patient in patient_list:
        try:
            insert_value = clinic_info[clinic_col][clinic_info['patient_id'] == patient]
            whole_info.loc[whole_info[whole_info['Tumor_Sample_Barcode'] == patient].index, clinic_col] = np.array(insert_value)
        except ValueError:
            whole_info.drop(whole_info[whole_info['Tumor_Sample_Barcode'] == patient].index, inplace=True)
            print (patient, ' is wrong.!!!!!!')
        num += 1
        print (num)
    print(whole_info)
    print ('Clinical infomation finished...')
    if os.path.exists(OLDER_PATH) == False:
        os.mkdir(OLDER_PATH)
    whole_info.to_csv(OLDER_PATH + 'whole.csv', index=False)
    for part in COMBINE_DICT:
        print(part, ' is collecting...')
        part_info = whole_info[whole_info['ICD_O3_site'].isin(COMBINE_DICT[part])]
        print(part_info)
        part_info.to_csv(OLDER_PATH + part +'.csv', index=False)
        print(part, ' has finished.')
