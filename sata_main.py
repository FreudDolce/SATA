# ========================================
# Used for satastic
# Written by Ji Hongchen
# 20200218
# ========================================

import pandas as pd
import numpy as np
import SATA_METHODS
import SATA_PRETREAT
import os
import shutil
import scipy
import argparse
import matplotlib.pyplot as plt
import datetime
import GETPERDATA
from pyfaidx import Fasta
from cfg import CFG


parser = argparse.ArgumentParser()
parser.add_argument('-f', help="the folder input and ouput.")
args = parser.parse_args()

cfg = CFG()

# The below items should be check, modify if nessary.
OUTPUT = cfg.output
CLINICAL_ITEM = cfg.clinical_item
MIN_ALIVE = cfg.min_alive
VOTE_CLASS = 'Class_4'
VOTE_LINE = 1
WITH_TYPICAL = 'no'
SATA_LIST = ['t_stage',
             'n_stage',
             'm_stage',
             'patient_gender',
             'patient_age',
             'patient_weight',
             'ajcc_stage',
             ]
# The above items should be check, modify if nessary.


now = datetime.datetime.now()
nowtime = now.strftime('%Y-%m-%d-%H-%M-%S')


if args.f:
    PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
        str(args.f) + '/analysis_data_' + str(args.f) + '/'
    if WITH_TYPICAL == 'no':
        PATH_LOG = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
            str(args.f) + '/analysis_log_' + str(args.f) + '/' + \
            str(nowtime) + '/'
    else:
        PATH_LOG = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
            str(args.f) + '/analysis_log_' + str(args.f) + '/' + \
            str(nowtime) + '-' + WITH_TYPICAL + '/'


def OutputFrame(dataframe, in_f, filename):
    output_path = PATH + in_f.split('.')[0] + '/'
    dataframe.to_csv(output_path + filename + '.csv')
    print('>>> ', filename, ' saved.')


def LoadFrame(ori_file, mode='b'):
    test_frame = pd.read_csv(ori_file, dtype={'patient_weight': str})
    test_frame = test_frame.astype('str')
    if mode == 'b':
        test_frame['Reference_Allele'] = test_frame['Reference_Allele'] + \
            test_frame['Tumor_Seq_Allele2']
        test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'GT'] = 'CA'
        test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'GC'] = 'CG'
        test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'GA'] = 'CT'
        test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'AT'] = 'TA'
        test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'TC'] = 'AG'
        test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'AC'] = 'TG'
        test_frame = SATA_PRETREAT.GetSurvTime(test_frame)
        test_frame = SATA_PRETREAT.NeuPreteate(test_frame)
        classed_frame = SATA_PRETREAT.Classifier(test_frame)
        classed_frame.drop(
            classed_frame[classed_frame['class_result'] == 0].index, inplace=True)
    classed_frame = SATA_PRETREAT.DropEmptyData(classed_frame,
                                                column='to_last_known_alive')
    return classed_frame


def CalSataItem(dataframe, outtxtfile, output_f, label='base'):
    sata_dict = cfg.clicfeat_dict
    for cal_item in SATA_LIST:
        calcframe = SATA_PRETREAT.DropEmptyData(dataframe, column=cal_item)
        print('{:*^120}'.format(cal_item + '_' + label), file=outtxtfile)
        print('{:*^120}'.format(cal_item + '_' + label))
        if sata_dict[cal_item] == 'ori':
            satasclass = SATA_METHODS.SataMethod(dataframe=calcframe,
                                                 class_col='class_result')
            try:
                chiasclass = satasclass.ChiSeqTest(chisq_col=cal_item)
                print(chiasclass, file=outtxtfile)
                OutputFrame(chiasclass, in_f=output_f,
                            filename=cal_item + '_' + label + '_as_class')
            except UnboundLocalError:
                print(cal_item, ' is loss', file=outtxtfile)
            sataasitem = SATA_METHODS.SataMethod(dataframe=calcframe,
                                                 class_col=cal_item)
            try:
                chiasitem = sataasitem.ChiSeqTest(chisq_col='class_result')
                print(chiasitem, file=outtxtfile)
                OutputFrame(chiasitem, in_f=output_f,
                            filename=cal_item + '_' + label + '_as_item')
            except UnboundLocalError:
                print(cal_item, ' is loss', file=outtxtfile)
        elif sata_dict[cal_item] == 'mid':
            sataasmid = SATA_METHODS.SataMethod(dataframe=calcframe)
            itemdescribe = sataasmid.Describe(descb_col=cal_item)
            itemttest = sataasmid.CalTTest(describeframe=itemdescribe)
            print(itemdescribe, file=outtxtfile)
            print(itemttest, file=outtxtfile)
            OutputFrame(itemdescribe, in_f=output_f,
                        filename=cal_item + '_' + label + '_describe')
            OutputFrame(itemttest, in_f=output_f,
                        filename=cal_item + '_' + label + '_ttest')
            if cal_item == 'to_last_known_alive':
                print('{:*^120}'.format('logrank'), file=outtxtfile)
                lr_frame = sataasmid.LogRankTest(
                    PATH + output_f.split('.')[0] + '/' + label + '_logrank')
                OutputFrame(lr_frame, in_f=output_f,
                            filename=label + '_logrank')
                print(lr_frame, file=outtxtfile)
                print(lr_frame)
        else:
            for i in sata_dict[cal_item]:
                for j in sata_dict[cal_item][i]:
                    calcframe[cal_item][calcframe[cal_item] == j] = i
            calcframe.drop(calcframe[calcframe[cal_item] == 'delete'].index,
                           inplace=True)
            satacaledcalss = SATA_METHODS.SataMethod(dataframe=calcframe,
                                                     class_col='class_result')
            try:
                chicaledclass = satacaledcalss.ChiSeqTest(chisq_col=cal_item)
                print(chicaledclass, file=outtxtfile)
                OutputFrame(chicaledclass, in_f=output_f,
                            filename=cal_item + '_' + label + '_as_class')
            except UnboundLocalError:
                print(cal_item, ' is loss', file=outtxtfile)
            satacaleditem = SATA_METHODS.SataMethod(dataframe=calcframe,
                                                    class_col=cal_item)
            try:
                chicaleditem = satacaleditem.ChiSeqTest(
                    chisq_col='class_result')
                print(chicaleditem, file=outtxtfile)
                OutputFrame(chicaleditem, in_f=output_f,
                            filename=cal_item + '_' + label + '_as_item')
            except UnboundLocalError:
                print(cal_item, ' is loss.', file=outtxtfile)
    print("{:=^120}".format('finish'))
    print("{:=^120}".format('finish'), file=outtxtfile)
    outtxtfile.close()


if __name__ == '__main__':
    if os.path.exists(PATH_LOG) == False:
        os.makedirs(PATH_LOG)
    pre_files = os.listdir(PATH)
    for i in pre_files:
        if '.csv' not in i:
            shutil.rmtree(PATH+i)
    files = os.listdir(PATH)
    for fi in files:
        output_path = PATH + fi.split('.')[0]
        if os.path.exists(output_path) == False:
            os.mkdir(output_path)
        overwritefile = open(PATH_LOG + fi.split('.')
                             [0] + '_base_ow.txt', 'a+')
        print('{:=^120}'.format(fi), file=overwritefile)
        print('{:=^120}'.format(fi))
        classedframe = LoadFrame(PATH + '/' + fi)
        if WITH_TYPICAL != 'no':
            classedframe = classedframe[classedframe['Hugo_Symbol']
                                        == WITH_TYPICAL]
        CalSataItem(dataframe=classedframe,
                    outtxtfile=overwritefile,
                    output_f=fi,
                    label='base')
        """
        if fi != 'whole.csv':
            voteframe = SATA_PRETREAT.ClassiVoter(dataframe=classedframe,
                                                  class_cal=VOTE_CLASS,
                                                  line=VOTE_LINE)
            voteframe.to_excel(
                output_path + '/voted_class.xlsx', index=False)
            try:
                wholeframe = pd.concat(
                    [wholeframe, voteframe], sort=True).reset_index(drop=True)
            except NameError:
                wholeframe = voteframe
            overwritefile = open(PATH_LOG + fi.split('.')
                                 [0] + '_vote_ow.txt', 'a+')
            print('{:=^120}'.format(fi), file=overwritefile)
            print('{:=^120}'.format(fi))
            print('>>>>>>>> VOTE WAY: ' + str(VOTE_CLASS) +
                  '>=' + str(VOTE_LINE), file=overwritefile)
            CalSataItem(dataframe=voteframe,
                        outtxtfile=overwritefile,
                        output_f=fi,
                        label='vote')
    if WITH_TYPICAL == 'no':
        overwritefile = open(PATH_LOG + 'whole_vote_ow.txt', 'a+')
        print('{:=^120}'.format('whole'), file=overwritefile)
        print('{:=^120}'.format('whole'))
        wholeframe.to_excel(PATH_LOG + 'whole_vote.xlsx', index=False)
        CalSataItem(dataframe=wholeframe,
                    outtxtfile=overwritefile,
                    output_f='whole.csv',
                    label='vote')
        """
