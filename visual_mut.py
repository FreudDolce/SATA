
# Used for visualise the mutation sequence
# Written by Ji Hongchen
# 20200419
# =========================================

from cfg import CFG
import argparse
import os
from pyfaidx import Fasta
import SATA_PRETREAT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib


cfg = CFG()
parser = argparse.ArgumentParser()
parser.add_argument('-f', help="the folder input and ouput.")
args = parser.parse_args()

GENEDICT = {}
VISUALDICT = {}
EXT = 10  # cfg.ext
CAMP = matplotlib.cm.hot
WHOLE_SEQUENCE = Fasta(
    '/Users/freud/Documents/MANU/lstmsom_data/clinical_data/hg38.fa')
CMAP = 'Reds'
NORM = matplotlib.colors.Normalize(vmin=0, vmax=1)
PATH = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
    str(args.f) + '/'
CALC_NUM = 3
CAL_NUM_CSV = 20


def VisiableSeq(test_frame):
    print('File_Loaded')
    test_frame = test_frame.astype('str')
    test_frame['Reference_Allele'] = test_frame['Reference_Allele'] + \
        test_frame['Tumor_Seq_Allele2']
    test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'GT'] = 'CA'
    test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'GC'] = 'CG'
    test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'GA'] = 'CT'
    test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'AT'] = 'TA'
    test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'AG'] = 'TC'
    test_frame['Reference_Allele'].loc[test_frame['Reference_Allele'] == 'AC'] = 'TG'
    print("********************************")
    # print(list(set(test_frame['Reference_Allele'])))
    test_frame = SATA_PRETREAT.GetSurvTime(test_frame)
    test_frame = SATA_PRETREAT.NeuPreteate(test_frame)
    classed_frame = SATA_PRETREAT.Classifier(test_frame)
    classed_frame.drop(
        classed_frame[classed_frame['class_result'] == 0].index, inplace=True)
    classed_frame.reset_index(drop=True, inplace=True)
    print('visiable frame loaded')
    pre_list = []
    post_list = []
    for i in range(len(classed_frame)):
        try:
            pre_list.append(str(WHOLE_SEQUENCE.get_seq(classed_frame.Chromosome[i],
                                                       int(
                                                           classed_frame.Start_Position[i])-EXT,
                                                       int(classed_frame.Start_Position[i])-1)).lower())
        except ValueError:
            print('!!!!!!!', i, ' wrong.')
            pre_list.append('nan')
        try:
            post_list.append(str(WHOLE_SEQUENCE.get_seq(classed_frame.Chromosome[i],
                                                        int(
                                                            classed_frame.Start_Position[i])+1,
                                                        int(classed_frame.Start_Position[i])+EXT)).lower())
        except ValueError:
            print('!!!!!!!', i, ' wrong.')
            post_list.append('nan')
        if i % 10000 == 0:
            print(i, ' has finished.')
    pre_list = np.array(pre_list)
    post_list = np.array(post_list)
    pre_list = pre_list.reshape(len(classed_frame), -1)
    post_list = post_list.reshape(len(classed_frame), -1)
    classed_frame['front_seq'] = pre_list
    classed_frame['behind_seq'] = post_list
    print(pd.crosstab(classed_frame['class_result'],
                      classed_frame['Reference_Allele']))
    return classed_frame


def MapVisiable(frame, has_whole='no'):
    colorpool = ['tomato', 'teal', 'gold', 'midnightblue']
    if has_whole != 'yes':
        frame['front_seq'] = frame['front_seq'] + \
            frame['Reference_Allele'] + frame['behind_seq']
    else:
        pass
    # print(frame['front_seq'])
    class_list = ['Class_1', 'Class_2', 'Class_3', 'Class_4', 'Class_5',
                  'Class_6', 'Class_7', 'Class_8', 'Class_9', 'Class_10']
    frame['calc_list'] = 'nothing'
    frame['calc_list'] = frame['front_seq'].apply(
        lambda x: x[20 - CALC_NUM: 22 + CALC_NUM]).tolist()
    for clas in class_list:
        numberframe = frame[frame['class_result'] == clas]
        total = numberframe['front_seq'].count()
        outframe = pd.DataFrame(columns=['count', 'dev'])
        outframe['count'] = numberframe['calc_list'].value_counts()[0: 20]
        outframe['dev'] = outframe['count'] / total
    fig = plt.figure(figsize=(2.5, 10), dpi=300)
    for index in class_list:
        dataframe = frame[frame['class_result'] == index]
        total_num = len(dataframe)
        nupseq = np.array(dataframe['front_seq'])
        nupseq = nupseq.reshape(-1, 1)
        wholeseq = np.zeros((len(dataframe), cfg.ext * 2 + 2)).astype(np.str)
        for i in range(len(wholeseq)):
            wholeseq[i] = list(nupseq[i][0])
        wholeseq = pd.DataFrame(wholeseq)
        vis_whole = np.zeros((4, 2*EXT+2))
        for i in range(2*EXT+2):
            try:
                a_upper = wholeseq[i].value_counts()['A'] / total_num
            except KeyError:
                a_upper = 0
            try:
                a_lower = wholeseq[i].value_counts()['a'] / total_num
            except KeyError:
                a_lower = 0
            vis_whole[0][i] = a_upper + a_lower
            try:
                t_upper = wholeseq[i].value_counts()['T'] / total_num
            except KeyError:
                t_upper = 0
            try:
                t_lower = wholeseq[i].value_counts()['t'] / total_num
            except KeyError:
                t_lower = 0
            vis_whole[1][i] = t_upper + t_lower
            try:
                c_upper = wholeseq[i].value_counts()['C'] / total_num
            except KeyError:
                c_upper = 0
            try:
                c_lower = wholeseq[i].value_counts()['c'] / total_num
            except KeyError:
                c_lower = 0
            vis_whole[2][i] = c_lower + c_upper
            try:
                g_upper = wholeseq[i].value_counts()['G'] / total_num
            except KeyError:
                g_upper = 0
            try:
                g_lower = wholeseq[i].value_counts()['g'] / total_num
            except KeyError:
                g_lower = 0
            vis_whole[3][i] = g_lower + g_upper
        drawframe = pd.DataFrame(vis_whole,
                                 index=['A', 'T', 'C', 'G'])
        ax = fig.add_subplot(10, 1, int(index.split('_')[1]))
        xLabel = np.arange(len(drawframe.columns)) + 1
        total_precent = [0] * len(drawframe.columns)
        for i in range(len(drawframe.index)):
            ax.bar(xLabel, np.array(drawframe.loc[drawframe.index[i]]),
                   bottom=total_precent,
                   label=str(drawframe.index[i]),
                   color=colorpool[i])
            total_precent = total_precent + \
                np.array(drawframe.loc[drawframe.index[i]])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels([])
        ax.set_xticks([-0.5, 23])
        #ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=10)
        #ax.legend(loc=(-1, 0))
        # ax.set_title(index)
        ax.set_xticklabels([])
    return (vis_whole, plt)


def CountMBforPatients(dataframe, calcitem='ICD_O3_site'):
    class_list = ['Class_1', 'Class_2', 'Class_3', 'Class_4', 'Class_5',
                  'Class_6', 'Class_7', 'Class_8', 'Class_9', 'Class_10']
    claclist = list(set(dataframe[calcitem]))
    print(claclist)
    countframe = pd.DataFrame(columns=class_list)
    for i in claclist:
        if dataframe[dataframe[calcitem] == i][calcitem].count() > 1000:
            countframe.loc[i] = 0
    for clas in class_list:
        classframe = dataframe[dataframe['class'] == clas]
        for c in claclist:
            if c in countframe.index:
                try:
                    countframe[clas].loc[c] = \
                        classframe[classframe[calcitem] ==
                                   c]['Tumor_Sample_Barcode'].value_counts().mean()
                except KeyError:
                    countframe[clas].loc[c] = 0
    return (countframe, countframe.index)


def RateMBforPatients(dataframe, ind, calcitem='ICD_O3_site'):
    class_list = ['Class_1', 'Class_2', 'Class_3', 'Class_4', 'Class_5',
                  'Class_6', 'Class_7', 'Class_8', 'Class_9', 'Class_10']
    rate_list = ['Rate' + class_list[i] for i in range(len(class_list))]
    rateframe = pd.DataFrame(columns=rate_list, index=ind)
    #for k in cfg.clicfeat_dict[calcitem]:
    #    dataframe[calcitem][dataframe[calcitem].isin(
    #        cfg.clicfeat_dict[calcitem][k])] = k
    for cr in rateframe.index:
        crframe = dataframe[dataframe[calcitem] == cr]
        for i in rate_list:
            rateframe[i].loc[cr] = crframe[i].mean()
    return rateframe


def MBScatter(countframe, rateframe):
    countframe = countframe.astype('int')
    countframe = countframe.apply(np.log2)
    fig = plt.figure(figsize=(4.8, 10), dpi=300)
    print(rateframe)
    for i in range(len(countframe.columns)):
        for j in range(len(countframe.index)):
            c_1 = 0.7 + 0.3 * np.log10(rateframe.iloc[j, i])
            if c_1 < 0:
                c_1 = 0
            c_3 = -0.35 * np.log10(rateframe.iloc[j, i])
            if c_3 > 1:
                c_3 = 1
            plt.scatter([i], [j],
                        s=11.5 * (countframe.iloc[j, i] + 1),
                        c=np.array([c_1, 0, c_3]).reshape(1, -1))
    plt.yticks(range(len(countframe.index)), countframe.index, fontsize=13)
    return plt


if __name__ == '__main__':
    fl = ['whole.csv']
    #fl = os.listdir(PATH + 'analysis_data_20200617/')
    #fl.remove('whole.csv')
    #countframe = np.array([1, 5, 10, 50, 100, 500])
    #countframe = pd.DataFrame(countframe)
    # print(countframe)
    #rateframe = np.array([0.005, 0.01, 0.05, 0.1, 0.5])
    #rateframe = pd.DataFrame(rateframe)
    # print(rateframe)
    for f in fl:
        if '.csv' in f:
            print(f.split('.')[0])
            countf = PATH + 'analysis_data_20200617/' + f
            ratef = PATH + 'analysis_patient_20200617/' + f
            wholeframe = pd.read_csv(countf)
            ratedf = pd.read_csv(ratef)
            # for p in cfg.clicfeat_dict['ICD_O3_pathology']:
            #    wholeframe['ICD_O3_pathology'].loc[wholeframe['ICD_O3_pathology'].isin(cfg.clicfeat_dict['ICD_O3_pathology'][p])] = p
            countframe, indexes = CountMBforPatients(
                wholeframe, calcitem='ICD_O3_site')
            countframe = countframe.sort_index(ascending=True)
            print(countframe)
            rateframe = RateMBforPatients(
                ratedf, ind=indexes, calcitem='ICD_O3_site')
            rateframe = rateframe.sort_index(ascending=True)
            print(rateframe)
            plt = MBScatter(countframe, rateframe)
            plt.tight_layout()
            plt.savefig(PATH + 'analysis_data_20200617/' +
                        f.split('.')[0] + '/site_scatter_1.tif')
    """
    file_list = os.listdir(PATH)
    file_list.remove('whole.csv')
    print(file_list)
    for f in file_list:
        print(f)
        if '.csv' in f:
            print('--------->', f, ' start...')
            visfile = PATH + f
            print('folder: ', visfile)
            visframe = VisiableSeq(visfile)
            partframe = VisiableSeq(visfile)
            try:
                totalvisframe = pd.concat([totalvisframe, partframe], ignore_index=True)
            except NameError:
                totalvisframe = partframe
            drawframe, plt= MapVisiable(frame=visframe)
            plt.tight_layout()
            plt.savefig(PATH + f.split('.')[0] + '/genepile.tif')
            print('--------->', f, 'finished')
    print(pd.crosstab(totalvisframe['class_result'], totalvisframe['Reference_Allele']))
    #totalvisframe['front_seq'] = totalvisframe['front_seq'] + \
    #                totalvisframe['Reference_Allele'] + totalvisframe['behind_seq']
    #class_list = ['Class_1', 'Class_2', 'Class_3', 'Class_4', 'Class_5',
    #                'Class_6', 'Class_7', 'Class_8', 'Class_9', 'Class_10']
    #for c in class_list:
    #    calcframe = totalvisframe[totalvisframe['class_result'] == c]
    #    outframe = calcframe['front_seq'].value_counts()
    #    outframe.to_csv(PATH + 'whole/' + c + '_base_count.csv')
    drawframe, plt = MapVisiable(frame=totalvisframe)
    plt.tight_layout()
    plt.savefig(PATH + 'whole/genepile.tif')
    """
