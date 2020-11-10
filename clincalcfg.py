# ===============================
# configure of patient satastic
# Written by Ji Hongchen
# 20200617
# ===============================

import pandas as pd
import numpy as np
import os


class CFG():
    def __init__(self):
        self.clicfeat_dict = {'ajcc_stage': {'I': ['Stage I', 'Stage IA', 'Stage IB', 'Stage 0'],
                                             'II': ['Stage II', 'Stage IIC', 'Stage IIB', 'Stage IIA'],
                                             'III': ['Stage IIIA',  'Stage III', 'Stage IIIC', 'Stage IIIB'],
                                             'IV': ['Stage IVC', 'Stage IV', 'Stage IVB', 'Stage IVA'],
                                             'delete': ['IS', 'ajcc_stage', 'data_empty', 'IhdII NOS', 'Stage X', 'data_loss']},
                              'ICD_O3_pathology': {'Epithelial carcinoma, NOS': [801, 802],
                                                   'Squamous cell carcinoma': [805, 807, 808],
                                                   'Basal cell carcinoma': [809],
                                                   'Transitional cell carcinoma': [812, 813],
                                                   'Adenocarcinoma': [814, 816, 817, 818, 820,
                                                                      821, 823, 824, 825, 826,
                                                                      829, 831, 833, 834, 835,
                                                                      837, 838],
                                                   'Adnexal and skin appendage carcinoma': [840],
                                                   'Cystic, mucinous and serous carcinoma': [844, 846,
                                                                                             848, 849],
                                                   'Ductal, lobular and medullary carcinoma': [850, 851,
                                                                                               852, 854],
                                                   'Acinar cell carcinoma': [855],
                                                   'Complex epithelial carcinoma': [856, 857],
                                                   'Thymic epithelial carcinoma': [858],
                                                   'Paragangliomas and glomus tumor': [868, 869, 870],
                                                   'Melanoma': [872, 873, 874, 877],
                                                   'Soft tissue tumor and sarcoma, NOS': [880],
                                                   'Fibromatous sarcoma': [881, 882, 883],
                                                   'Lipomatous sarcoma': [885],
                                                   'Myomatous sarcoma': [889],
                                                   'Complex mixed and stromal sarcoma': [895, 898],
                                                   'Phyllodes tumor, malignant': [902],
                                                   'Synovioma sarcoma': [904],
                                                   'Mesothilioma, malignant': [905],
                                                   'Germ cell tumor, malignant': [906, 907, 908],
                                                   'Gliomas': [938, 940, 944, 945],
                                                   'Nerve sheath tumors': [954],
                                                   'Malignant lymphoma': [968]},
                              't_stage': {'T1': ['T1b', 'T1a1', 'T1', 'T1b1', 'T1a', 'T1c', 'T1b2'],
                                          'T2': ['T2', 'T2a', 'T2a2', 'T2b', 'T2c', 'T2a1', ],
                                          'T3': ['T3b', 'T3a', 'T3', 'T3c'],
                                          'T4': ['T4d', 'T4c', 'T4', 'T4b', 'T4a', 'T4e'],
                                          'delete': ['data_empty', 'data_loss', 'TX', 'T0', 'Tis']},
                              'n_stage': {'N0': ['N0', 'N0 (mol+)', 'N0 (i+)', 'N0 (i-)'],
                                          'N1': ['N1mi', 'N1a', 'N1', 'N1b', 'N1c'],
                                          'N2': ['N2a', 'N2', 'N2c', 'N2b'],
                                          'N3': ['N3a', 'N3b', 'N3', 'N3c'],
                                          'delete': ['data_empty', 'empty_data', 'NX']},
                              'm_stage': {'M0': ['cM0 (i+)', 'M0'],
                                          'M1': ['M1a', 'M1', 'M1b', 'M1c'],
                                          'delete': ['MX', 'data_loss', 'data_empty']},
                              'patient_gender': {'Male': ['MALE'], 'Female': ['FEMALE'],
                                                 'delete': ['data_loss', 'data_empty']}}

        self.clicfeat_dict['patient_age'] = {}
        self.clicfeat_dict['patient_weight'] = {}
        for i in range(2, 18):
            self.clicfeat_dict['patient_age'][str(i * 5 + 1) + ' ~ ' + str(i * 5 + 5)] =\
                [str(j) for j in range(i * 5 + 1, i * 5 + 6)]
        self.clicfeat_dict['patient_age']['delete'] = [
            'data_loss', 'data_empty']
        self.clicfeat_dict['patient_weight']['0 ~ 40'] = [
            str(l) for l in range(0, 41)]
        for k in range(4, 12):
            self.clicfeat_dict['patient_weight'][str(k * 10 + 1) + ' ~ ' + str(k * 10 + 10)] =\
                [str(h) for h in range(k * 10 + 1, k * 10 + 11)]
        self.clicfeat_dict['patient_weight']['>120'] = [
            str(a) for a in range(120, 300)]
        self.clicfeat_dict['patient_weight']['delete'] = [
            'data_loss', 'data_empty']
