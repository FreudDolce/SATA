# =========================================================
# Used for get mutation data and sequences data
# Written by Ji Hongchen
# 20191101
# =========================================================

import numpy as np
import pandas as pd
from pyfaidx import Fasta
import os

WHOLE_SEQUENCE = '/media/ji/data/whole_sequence/hg38.fa'
TEST_MAF = '/media/ji/data/programme/GENE_ANALYSIS/examples/TCGA.CHOL.somaticsniper.5d5dfadb-f18d-4f19-aff4-166dac7b92df.DR-10.0.somatic.maf'
TARGET_FOLDER = r'/media/ji/data/mut_info/'


def GetWholeSequence(seq_path):
    """
    Return the whole
    """
    genes = Fasta(seq_path)
    return genes


def MutDict():
    """
    return the code dict of genes
    """
    mutdict = {
        'A': [1.0],
        'T': [2.0],
        'C': [3.0],
        'G': [4.0],
        'N': [0.0]
    }
    return mutdict


def _get_patient_id(name):
    p_name = name.split('-')
    p_id = '-'.join(p_name[0:3])
    return p_id


def ReadMafData(maf_file):
    kwargs = {'sep': '\t', 'comment': '#'}
    mut_info = pd.read_csv(maf_file, **kwargs)
    mut_list = []
    for i in range(len(mut_info)):
        mut_list.append([_get_patient_id(mut_info['Tumor_Sample_Barcode'][i]),
                         mut_info['Hugo_Symbol'][i],
                         mut_info['Chromosome'][i],
                         mut_info['Start_Position'][i],
                         mut_info['End_Position'][i],
                         mut_info['Reference_Allele'][i],
                         mut_info['Tumor_Seq_Allele1'][i],
                         mut_info['Tumor_Seq_Allele2'][i]])
    mut_list = np.array(mut_list)
    return mut_list


def GetMutSeq(geneseq, mutpos, extend=25):
    """
    parameters:geneseq, mutpos
    geneseq:genes
    mutpos:[chr, start_position, ref_all, t_all_1, t_all_2]
    return (pre_list(extend before mut), post_list(extend after mut))
    """
    mutdict = MutDict()
    pre_list = str(geneseq.get_seq(
        mutpos[0], mutpos[1]-extend, mutpos[1]-1)).upper()
    post_list = str(geneseq.get_seq(
        mutpos[0], mutpos[1]+1, mutpos[1]+extend)).upper()
    post_list = post_list[:: -1]
    pre_list_mat = np.zeros((extend+1, 1))
    post_list_mat = np.zeros((extend+1, 1))
    mut_mat = np.array([mutdict[mutpos[2][0]],
                        mutdict[mutpos[3][0]],
                        mutdict[mutpos[4][0]]]).reshape(1, -1)

    for base in pre_list:
        gene_mat = np.array(mutdict[base]).reshape(1, -1)
        gene_mat = np.hstack((gene_mat, gene_mat, gene_mat))
        if pre_list_mat.sum() == 0:
            pre_list_mat = gene_mat
        else:
            pre_list_mat = np.vstack((pre_list_mat, gene_mat))
    pre_list_mat = np.r_[pre_list_mat, mut_mat]

    for post_base in post_list:
        gene_mat = np.array(mutdict[post_base]).reshape(1, -1)
        gene_mat = np.hstack((gene_mat, gene_mat, gene_mat))
        if post_list_mat.sum() == 0:
            post_list_mat = gene_mat
        else:
            post_list_mat = np.vstack((post_list_mat, gene_mat))
    post_list_mat = np.r_[post_list_mat, mut_mat]

    return (pre_list_mat, post_list_mat)


if __name__ == '__main__':
    folders = os.listdir('/media/ji/data/TCGA_DATA/MAF/')
    outputtlist = None
    for folder in folders:
        files = os.listdir('/media/ji/data/TCGA_DATA/MAF/' + folder)
        for file in files:
            if ('somaticsniper' in file):
                mutlist = ReadMafData(
                    '/media/ji/data/TCGA_DATA/MAF/' + folder + '/' + file)
                print(file + 'listed.')
                print(
                    '===========================================================================')
                if outputtlist is None:
                    outputtlist = mutlist
                else:
                    outputtlist = np.vstack((outputtlist, mutlist))
                print(outputtlist.shape)
    len_outputlist = len(outputtlist)
    col_name = ['Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position',
                'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
    for i in range(len_outputlist):
        print(outputtlist[i])
        lucky_num = np.random.choice(len_outputlist)
        file_num = int(lucky_num // 5000)  # 5000 items per file
        print('lucky_num=' + str(file_num))
        try:
            existed_list = pd.read_csv(
                TARGET_FOLDER + 'mut_info_' + str(file_num) + '.csv')
            existed_list = np.array(existed_list)
            output_list = np.vstack((existed_list, outputtlist[i]))
            output_list = pd.DataFrame(output_list)
            output_list.to_csv(TARGET_FOLDER + 'mut_info_' + str(file_num) + '.csv',
                               header=col_name,
                               index=None)
        except FileNotFoundError:
            output_list = pd.DataFrame([outputtlist[i]])
            output_list.to_csv(TARGET_FOLDER + 'mut_info_' + str(file_num) + '.csv',
                               header=col_name,
                               index=None)
