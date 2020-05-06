# =========================================
# Used for visualise the mutation sequence
# Written by Ji Hongchen
# 20200419
# =========================================

import os
from pyfaidx import Fasta
import SATA_PRETREAT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

GENEDICT = {}
VISUALDICT = {}
VISIABLEFOLD = '20200305'
ACCORDCOL = 'pred_1'
THRESH = [3, 1.5, 1, 0]
EXT = 15
CAMP = matplotlib.cm.hot
WHOLE_SEQUENCE = Fasta(
    '/Users/freud/Documents/MANU/lstmsom_data/clinical_data/hg38.fa')
CMAP = 'OrRd'
NORM = matplotlib.colors.Normalize(vmin=0, vmax=1)


def VisiableSeq(visiablefile):
    test_frame = pd.read_csv(visiablefile)
    print('File_Loaded')
    COL_NAME = ['gene_name', 'mut_pos', 'mut_gene',
                'Chromosome', 'class_name', 'front_seq', 'behind_seq']
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
    classed_frame.reset_index(drop=True, inplace=True)
    zerosframe = np.zeros((len(classed_frame), len(COL_NAME)))
    visiableframe = pd.DataFrame(zerosframe, columns=COL_NAME)
    visiableframe.gene_name = classed_frame.Hugo_Symbol
    visiableframe.Chromosome = classed_frame.Chromosome
    visiableframe.mut_pos = classed_frame.Start_Position
    visiableframe.mut_gene = classed_frame.Reference_Allele
    visiableframe.class_name = classed_frame.class_result
    print('visiable frame loaded')
    for i in range(len(visiableframe)):
        try:
            pre_list = str(WHOLE_SEQUENCE.get_seq(visiableframe.Chromosome[i],
                                                  int(visiableframe.mut_pos[i])-EXT,
                                                  int(visiableframe.mut_pos[i])-1)).upper()
            post_list = str(WHOLE_SEQUENCE.get_seq(visiableframe.Chromosome[i],
                                                   int(visiableframe.mut_pos[i])+1,
                                                   int(visiableframe.mut_pos[i])+EXT)).upper()
        except ValueError:
            print(i)
            print(visiableframe.loc[i])
        if i % 10000 == 0:
            print(i, ' has finished.')
        visiableframe.front_seq[i] = pre_list
        visiableframe.behind_seq[i] = post_list
    print(visiableframe)
    return visiableframe


def MapVisiable(frame, output_path):
    frame['front_seq'] = frame['front_seq'] + \
        frame['mut_gene'] + frame['behind_seq']
    class_list = list(set(np.array(frame['class_name'])))
    for index in class_list:
        dataframe = frame[frame['class_name'] == index]
        dataframe.reset_index(drop=True, inplace=True)
        total_num = len(dataframe)
        wholeseq = ''.join(list(np.array(dataframe['front_seq'])))
        wholeseq = list(wholeseq)
        wholeseq = np.array(wholeseq)
        wholeseq = wholeseq.reshape(-1, 2*EXT+2)
        wholeseq = pd.DataFrame(wholeseq)
        vis_whole = np.zeros((2*EXT+2, 4))
        for i in range(2*EXT+2):
            try:
                vis_whole[i][0] = wholeseq[i].value_counts()['A'] / total_num
            except KeyError:
                vis_whole[i][0] = 0
            try:
                vis_whole[i][1] = wholeseq[i].value_counts()['T'] / total_num
            except KeyError:
                vis_whole[i][1] = 0
            try:
                vis_whole[i][2] = wholeseq[i].value_counts()['C'] / total_num
            except KeyError:
                vis_whole[i][2] = 0
            try:
                vis_whole[i][3] = wholeseq[i].value_counts()['G'] / total_num
            except KeyError:
                vis_whole[i][3] = 0
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(vis_whole, cmap=CMAP, norm=NORM)
        plt.colorbar(im)
        plt.savefig(output_path + '/hot_' + str(index) + '.png')
        plt.close()


if __name__ == '__main__':
    visiablefolder = '/Users/freud/Documents/MANU/lstmsom_data/exp' + \
        VISIABLEFOLD + '/analysis_data_' + VISIABLEFOLD + '/'
    file_list = os.listdir(visiablefolder)
    file_list.remove('whole.csv')
    print(file_list)
    for f in file_list:
        print(f)
        if '.csv' in f:
            print('--------->', f, ' start...')
            visfile = visiablefolder + f
            out_file = visiablefolder + f.split('.')[0] + '/visable.csv'
            print('folder: ', visfile)
            print('Output in: ', out_file)
            visframe = VisiableSeq(visfile)
            MapVisiable(visframe, visiablefolder+f.split('.')[0])
            visframe.to_csv(out_file, index=0)
            print('--------->', f, 'finished')
