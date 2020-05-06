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


parser = argparse.ArgumentParser()
parser.add_argument('-f', help="the folder input and ouput.")
args = parser.parse_args()


# The below items should be check, modify if nessary.
ACCORDCOL = 'pred_1'    # class according to which column
THRESH = [3, 1.5, 1, 0]            # boundary of classes
OUTPUT = True
MODE = 'p_b'
# 'p_b: base and patient
# 'pre': return a dict as {class:number}
# 'b': base only
# 'p': patient only, shoule give parameter "mode_p_frame": the clinical frame
VOTE_LINE = 0.1                     # The ratio to decide clinical class
# if vote_method = 'fold', use FREQ_BIAS[0]
FREQ_BIAS = [0, 0]
VOTE_METHOD = 'fold'                # quadrant or fold
# Use when mode "quadrant"
VOTE_DIRECT = np.array([[0, 1], [1, 0], [0, -1], [-1, 0]])
CLINICAL_ITEM = ['Tumor_Sample_Barcode', 'ICD_O3_pathology', 'ICD_O3_site',
                 'ajcc_stage', 'has_family_history', 'has_risk_factor',
                 'has_vasinvation', 'last_tumor_statu', 'm_stage', 'n_stage',
                 'neo_adj_info', 'patient_age', 'patient_gender', 'patient_height',
                 'patient_race', 'patient_status', 'patient_weight', 'to_last_known_alive',
                 'vote_class']
DRAW_MARK = [1, 1, -1, -1]
MIN_ALIVE = 15
# The above items should be check, modify if nessary.

if args.f:
    PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
        str(args.f) + '/analysis_data_' + str(args.f) + '/'


def OutputFrame(dataframe, in_f, filename):
    output_path = PATH + in_f.split('.')[0] + '/'
    dataframe.to_csv(output_path + filename + '.csv')
    print('>>> ', filename, ' saved.')


def AdjustFreqDict(frequence_dict, bias=FREQ_BIAS):
    adjusted_dict = {}
    total = 0
    for item in list(frequence_dict):
        total += frequence_dict[item]
    mean = total/len(frequence_dict)
    for j in list(frequence_dict):
        adjusted_dict[j] = frequence_dict[j] - \
            (frequence_dict[j] - mean) * frequence_dict[j] * bias[0]
    return adjusted_dict


def DrawPatientClass(dataframe, draw_mark, path):
    if draw_mark[3] == 0:
        dataframe['add'] = 20
    drawframe = np.array(dataframe)
    _y = draw_mark[0] * drawframe[:, 0] + draw_mark[2] * drawframe[:, 2]
    _x = draw_mark[1] * drawframe[:, 1] + draw_mark[3] * drawframe[:, 3]
    plt.scatter(_x, _y,
                c='r', s=5,
                alpha=0.2)
    plt.savefig(path+'/patient_class.png')
    plt.close()


def CalResult(f, fq_dict={}, mode=MODE, HAS_OUTPUT=OUTPUT, mode_p_frame=np.array([])):
    output_path = PATH + f.split('.')[0]
    if os.path.exists(output_path) == False:
        os.mkdir(output_path)
    now = datetime.datetime.now()
    nowtime = now.strftime('%Y-%m-%d-%H-%M-%S')
    overview_file = open(output_path+'/'+'overview_'+str(nowtime)+'.txt', 'a+')
    print('======================================================================',
          file=overview_file)
    print('Run time: ', now, file=overview_file)
    print('Using experience: ', args.f, file=overview_file)
    print('>>>'+str(f), file=overview_file)
    print('>>>vote method: ', VOTE_METHOD, file=overview_file)
    print('Vote direct: ', VOTE_DIRECT, file=overview_file)
    print('>>>vote thresh: ', VOTE_LINE, file=overview_file)
    print('>>>vote bias: ', FREQ_BIAS, file=overview_file)
    print('>>>min alive calculated: ', MIN_ALIVE, file=overview_file)
    if mode != 'p':
        test_frame = pd.read_csv(PATH+f)
        test_frame = test_frame.astype('str')
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
        classed_frame = SATA_PRETREAT.Classifier(test_frame,
                                                 accordcol=ACCORDCOL,
                                                 thresh=THRESH)
        classed_frame.drop(
            classed_frame[classed_frame['class_result'] == 0].index, inplace=True)
        classed_frame = SATA_PRETREAT.DropEmptyData(classed_frame,
                                                    column='to_last_known_alive')
        sata_methods = SATA_METHODS.SataMethod(dataframe=classed_frame,
                                               class_col='class_result')
        print('>>>total patient_number: '+str(len(np.unique(np.array(classed_frame['Tumor_Sample_Barcode'])))),
              file=overview_file)
        print('======================================================================',
              file=overview_file)
        print('************************-> Mutation mode  <-**************************',
              file=overview_file)
        chi2_mutmode = sata_methods.ChiSeqTest(chisq_col='Reference_Allele')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_mutmode, f, 'chi2_mutmode')
        print(chi2_mutmode, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************->  allele1 base  <-**************************',
              file=overview_file)
        chi2_allele1 = sata_methods.ChiSeqTest(chisq_col='Tumor_Seq_Allele1')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_allele1, f, 'chi2_allele1')
        print(chi2_allele1, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************->  allele2 base  <-**************************',
              file=overview_file)
        chi2_allele2 = sata_methods.ChiSeqTest(chisq_col='Tumor_Seq_Allele2')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_allele2, f, 'chi2_allele2')
        print(chi2_allele2, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('*********************-> Chromosome distribute  <-*********************',
              file=overview_file)
        chi2_chr = sata_methods.ChiSeqTest(chisq_col='Chromosome')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_chr, f, 'chi2_chr')
        print(chi2_chr, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************ -> Patient age  < -*************************',
              file=overview_file)
        age_descrb = sata_methods.Describe(descb_col='patient_age')
        age_ttest = sata_methods.CalTTest(describeframe=age_descrb)
        if HAS_OUTPUT == True:
            OutputFrame(age_descrb, f, 'age_descrb')
            OutputFrame(age_ttest, f, 'age_ttest')
        print(age_descrb, file=overview_file)
        print(age_ttest, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************-> Patient weight  <-*************************',
              file=overview_file)
        weight_descrb = sata_methods.Describe(descb_col='patient_weight')
        weight_ttest = sata_methods.CalTTest(describeframe=weight_descrb)
        if HAS_OUTPUT == True:
            OutputFrame(weight_descrb, f, 'weight_descrb')
            OutputFrame(weight_ttest, f, 'weight_ttest')
        print(weight_descrb, file=overview_file)
        print(weight_ttest, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************-> Patient gender  <-*************************',
              file=overview_file)
        chi2_gender = sata_methods.ChiSeqTest(chisq_col='patient_gender')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_gender, f, 'chi2_gender')
        print(chi2_gender, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************->  Patient race   <-*************************',
              file=overview_file)
        chi2_race = sata_methods.ChiSeqTest(chisq_col='patient_race')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_race, f, 'chi2_race')
        print(chi2_race, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************-> ICD O3 Phology <-**************************',
              file=overview_file)
        chi2_icdpathology = sata_methods.ChiSeqTest(
            chisq_col='ICD_O3_pathology')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_icdpathology, f, 'chi2_icdpathology')
        print(chi2_icdpathology, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************->  ICD O3 site   <-**************************',
              file=overview_file)
        chi2_icdsite = sata_methods.ChiSeqTest(chisq_col='ICD_O3_site')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_icdsite, f, 'chi2_icdsite')
        print(chi2_icdsite, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************-> AJCC stage 2017 <-*************************',
              file=overview_file)
        chi2_ajccstage = sata_methods.ChiSeqTest(chisq_col='ajcc_stage')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_ajccstage, f, 'chi2_ajccstage')
        print(chi2_ajccstage, file=overview_file)
        print('***********  T stage', file=overview_file)
        chi2_tstage = sata_methods.ChiSeqTest(chisq_col='t_stage')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_tstage, f, 'chi2_tstage')
        print(chi2_tstage, file=overview_file)
        print('***********  N stage', file=overview_file)
        chi2_nstage = sata_methods.ChiSeqTest(chisq_col='n_stage')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_nstage, f, 'chi2_nstage')
        print(chi2_nstage, file=overview_file)
        print('***********  M stage', file=overview_file)
        chi2_mstage = sata_methods.ChiSeqTest(chisq_col='m_stage')
        if HAS_OUTPUT == True:
            OutputFrame(chi2_mstage, f, 'chi2_mstage')
        print(chi2_mstage, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************-> alive describe <-**************************',
              file=overview_file)
        descrb_alive = sata_methods.Describe(descb_col='to_last_known_alive')
        print(descrb_alive, file=overview_file)
        ttest_alive = sata_methods.CalTTest(describeframe=descrb_alive)
        if HAS_OUTPUT == True:
            OutputFrame(descrb_alive, f, 'alive_describe')
            OutputFrame(ttest_alive, f, 'ttest_alive')
        print(ttest_alive, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('************************-> logrank result <-**************************',
              file=overview_file)
        lr_frame = sata_methods.LogRankTest(
            outpath=PATH+f.split('.')[0]+'/'+'logrank')
        if HAS_OUTPUT == True:
            OutputFrame(lr_frame, f, 'logistic_regression')
        print(lr_frame, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
    if mode == 'pre':
        _frequence_dict = {}
        total_bases = descrb_alive.loc['count'].sum()
        for c in descrb_alive.columns:
            _frequence_dict[c] = descrb_alive[c]['count'] / total_bases
        print('----------------------------------------------------------------------',
              file=overview_file)
        print('frequence_dict: ', _frequence_dict, file=overview_file)
        output_dict = AdjustFreqDict(frequence_dict=_frequence_dict,
                                     bias=FREQ_BIAS)
        print('adjusted_dict: ', output_dict, file=overview_file)
        return output_dict
    elif ((mode == 'p_b') or (mode == 'p')):
        print('----------------------------------------------------------------------',
              file=overview_file)
        print('! Below is the calculate used clinical class !', file=overview_file)
        print('----------------------------------------------------------------------',
              file=overview_file)
        if len(mode_p_frame) == 0:
            vote_frame, calc_frame = SATA_PRETREAT.ClassiVoter(dataframe=classed_frame,
                                                               freq_dict=fq_dict,
                                                               method=VOTE_METHOD,
                                                               vote_line=VOTE_LINE,
                                                               direct=VOTE_DIRECT,
                                                               bias=FREQ_BIAS)
        else:
            vote_frame = mode_p_frame[0]
            calc_frame = mode_p_frame[1]
        vote_frame.drop(
            vote_frame[vote_frame['to_last_known_alive'] <= MIN_ALIVE].index, inplace=True)
        DrawPatientClass(calc_frame,
                         draw_mark=DRAW_MARK,
                         path=output_path)
        output_vote_frame = vote_frame[CLINICAL_ITEM]
        output_vote_frame.to_excel(output_path+'/'+'clinical_data.xlsx')
        print('>>> Analysis patient: ', len(
            vote_frame), '.', file=overview_file)
        calc_frame.to_csv(output_path+'/'+'patient_class.csv')
        sata_method = SATA_METHODS.SataMethod(dataframe=vote_frame,
                                              class_col='vote_class')
        print('****************** -> Voted Patient age  < -**********************',
              file=overview_file)
        vote_age_descrb = sata_method.Describe(descb_col='patient_age')
        vote_age_ttest = sata_method.CalTTest(describeframe=vote_age_descrb)
        if HAS_OUTPUT == True:
            OutputFrame(vote_age_descrb, f, 'vote_age_descrb')
            OutputFrame(vote_age_ttest, f, 'vote_age_ttest')
        print(vote_age_descrb, file=overview_file)
        print(vote_age_ttest, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('**********************-> Vote Patient weight  <-**********************',
              file=overview_file)
        vote_weight_descrb = sata_method.Describe(descb_col='patient_weight')
        vote_weight_ttest = sata_method.CalTTest(
            describeframe=vote_weight_descrb)
        if HAS_OUTPUT == True:
            OutputFrame(vote_weight_descrb, f, 'vote_weight_descrb')
            OutputFrame(vote_weight_ttest, f, 'vote_weight_ttest')
        print(vote_weight_descrb, file=overview_file)
        print(vote_weight_ttest, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('**********************-> Vote Patient gender  <-**********************',
              file=overview_file)
        vote_chi2_gender = sata_method.ChiSeqTest(chisq_col='patient_gender')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_gender, f, 'vote_chi2_gender')
        print(vote_chi2_gender, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('***********************-> Vote Patient race <-************************',
              file=overview_file)
        vote_chi2_race = sata_method.ChiSeqTest(chisq_col='patient_race')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_race, f, 'vote_chi2_race')
        print(vote_chi2_race, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('**********************-> Vote ICD O3 Phology <-***********************',
              file=overview_file)
        vote_chi2_icdpathology = sata_method.ChiSeqTest(
            chisq_col='ICD_O3_pathology')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_icdpathology, f, 'vote_chi2_icdpathology')
        print(vote_chi2_icdpathology, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('**********************-> Vote ICD O3 site  <-*************************',
              file=overview_file)
        vote_chi2_icdsite = sata_method.ChiSeqTest(chisq_col='ICD_O3_site')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_icdsite, f, 'vote_chi2_icdsite')
        print(vote_chi2_icdsite, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('*********************-> Vote AJCC stage 2017 <-***********************',
              file=overview_file)
        vote_chi2_ajccstage = sata_method.ChiSeqTest(chisq_col='ajcc_stage')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_ajccstage, f, 'vote_chi2_ajccstage')
        print(vote_chi2_ajccstage, file=overview_file)
        print('***********  T stage', file=overview_file)
        vote_chi2_tstage = sata_method.ChiSeqTest(chisq_col='t_stage')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_tstage, f, 'vote_chi2_tstage')
        print(vote_chi2_tstage, file=overview_file)
        print('***********  N stage', file=overview_file)
        vote_chi2_nstage = sata_method.ChiSeqTest(chisq_col='n_stage')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_nstage, f, 'vote_chi2_nstage')
        print(vote_chi2_nstage, file=overview_file)
        print('***********  M stage', file=overview_file)
        vote_chi2_mstage = sata_method.ChiSeqTest(chisq_col='m_stage')
        if HAS_OUTPUT == True:
            OutputFrame(vote_chi2_mstage, f, 'vote_chi2_mstage')
        print(vote_chi2_mstage, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('**********************-> Vote alive describe <-***********************',
              file=overview_file)
        vote_descrb_alive = sata_method.Describe(
            descb_col='to_last_known_alive')
        print(vote_descrb_alive, file=overview_file)
        vote_ttest_alive = sata_method.CalTTest(
            describeframe=vote_descrb_alive)
        if HAS_OUTPUT == True:
            OutputFrame(vote_descrb_alive, f, 'vote_alive_describe')
            OutputFrame(vote_ttest_alive, f, 'vote_ttest_alive')
        print(vote_ttest_alive, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        print('**********************-> Vote logrank result <-***********************',
              file=overview_file)
        vote_lr_frame = sata_method.LogRankTest(
            outpath=PATH+f.split('.')[0]+'/'+'vote_logrank')
        if HAS_OUTPUT == True:
            OutputFrame(vote_lr_frame, f, 'vote_logistic_regression')
        print(vote_lr_frame, file=overview_file)
        print('**********************************************************************',
              file=overview_file)
        return (vote_frame, calc_frame)
    overview_file.close()
    print('======================================================================')
    print('***********************', f, ' finshed ************************')
    print('======================================================================')


if __name__ == '__main__':
    pre_files = os.listdir(PATH)
    for i in pre_files:
        if '.csv' not in i:
            shutil.rmtree(PATH+i)
    files = os.listdir(PATH)
    files.remove('whole.csv')
    print('======================================================================')
    print('>>> Calculating total bases:')
    print('======================================================================')
    adjust_dict = CalResult('whole.csv', mode='pre')
    for fi in files:
        print(fi)
        print('======================================================================')
        print('*************************** ',
              fi, ' ***************************')
        print('======================================================================')
        voteframe, calcframe = CalResult(fi,
                                         fq_dict=adjust_dict,
                                         mode='p_b')
        try:
            wholeframe = pd.concat([wholeframe, voteframe])
            wholecalcframe = pd.concat([wholecalcframe, calcframe])
        except NameError:
            wholeframe = voteframe
            wholecalcframe = calcframe
    print('======================================================================')
    print('******************************* whole  *******************************')
    print('======================================================================')
    CalResult('whole.csv',
              fq_dict=adjust_dict,
              mode='p',
              mode_p_frame=(wholeframe, wholecalcframe))
