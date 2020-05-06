# ==================================================
# Used for combine files according to ICD code
# Written by Ji Hongchen
# 20200309
# ==================================================

import pandas as pd
import os
import shutil
import argparse

COMBINE_DICT = {'head_neck':    ['00', '01', '02', '03', '04',
                                 '05', '06', '07', '08', '09',
                                 '10', '13', '14'],
                'esophagus':    ['15'],
                'stomach_si':   ['16', '17'],
                'colon':        ['18'],
                'rectom_anus':  ['19', '20', '21'],
                'liver':        ['22'],
                'baliary':      ['24'],
                'pancerase':    ['25'],
                'nasal':        ['30', '32'],
                'lung_bron':    ['34'],
                'thymus':       ['37'],
                'heart_pleura': ['38'],
                'bone_joints':  ['40', '41'],
                'skin':         ['44'],
                'peritoneum':   ['48'],
                'breast':       ['50'],
                'girl_genital': ['51', '52', '53'],
                'uters_uteri':  ['54', '55'],
                'ovary':        ['56'],
                'prost_gland':  ['61'],
                'testis':       ['62'],
                'kidney':       ['64'],
                'bladder':      ['67'],
                'eyes_adnexa':  ['69'],
                'meni_brain':   ['70', '71'],
                'thyroid':      ['73'],
                'adrenal':      ['74'],
                'lymph_node':   ['77']}

parser = argparse.ArgumentParser()
parser.add_argument('-e', help='the folder need to treated.')
args = parser.parse_args()

FOLDER = '/Users/freud/Documents/MANU/lstmsom_data/'

if args.e:
    EXP_NUM = args.e
    PATH = FOLDER + 'exp' + EXP_NUM + '/'

if __name__ == '__main__':
    os.makedirs(PATH + 'analysis_data_' + EXP_NUM)
    for part in COMBINE_DICT:
        serial_list = COMBINE_DICT[part]
        if len(serial_list) == 1:
            shutil.copyfile(PATH + 'cluster_result_' + EXP_NUM + '/flag_C' + serial_list[0] + '.csv',
                            PATH + 'analysis_data_' + EXP_NUM + '/' + part + '.csv')
            print(part, ' combinded.')
        else:
            partframe_1 = pd.read_csv(PATH + 'cluster_result_' + EXP_NUM + '/flag_C' + serial_list[0] + '.csv')
            for i in range(1, len(serial_list)):
                partframe_2 = pd.read_csv(PATH + 'cluster_result_' + EXP_NUM + '/flag_C' + serial_list[i] + '.csv')
                partframe_1 = pd.concat([partframe_1, partframe_2])
            print(part, ' combinded')
            partframe_1.to_csv(PATH + 'analysis_data_' + EXP_NUM + '/' + part + '.csv', index=False)
    print ('--------------------> finished <---------------------')
