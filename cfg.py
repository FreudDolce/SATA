# ========================================
# Configuration of SATA
# Written by Ji Hongchen
# 20200508
# ========================================

import numpy as np


class CFG():

    def __init__(self):
        # self.accordcol = 'pred_1'
        # self.thresh = [2, 0.7, 0.05, -0.2, -0.5, -0.8, -1.75]
        self.output = True
        self.mode_as_class = 'p_b'
        self.clinical_item = ['Tumor_Sample_Barcode', 'ICD_O3_pathology', 'ICD_O3_site',
                              'ajcc_stage',  'has_risk_factor', 'Chromosome',
                              'has_vasinvation', 'last_tumor_statu', 'm_stage', 'n_stage',
                              'patient_age', 'patient_gender', 'patient_height',
                              'patient_race', 'patient_status', 'patient_weight', 'to_last_known_alive', 'class_result']
        self.min_alive = 15
        self.mode_as_clicfeature = 'second_class'
        self.calc_clic_feature = ['ICD_O3_pathology', 'ICD_O3_site', 'ajcc_stage',
                                  'patient_age', 'patient_gender', 'patient_race',
                                  'patient_weight', 'to_last_known_alive', 't_stage', 'n_stage', 'm_stage']
        self.clicfeat_dict = {'ICD_O3_pathology': {'Epithelial carcinoma, NOS': [801, 802],
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
                              'ICD_O3_site': 'ori',
                              'Chromosome': 'ori',
                              'ajcc_stage': {'III': ['Stage IIIA',  'Stage III', 'Stage IIIC', 'Stage IIIB'],
                                             'I': ['Stage I', 'Stage IA', 'Stage IB', 'Stage 0'],
                                             'II': ['Stage II', 'Stage IIC', 'Stage IIB', 'Stage IIA'],
                                             'IV': ['Stage IVC', 'Stage IV', 'Stage IVB', 'Stage IVA'],
                                             'delete': ['IS', 'ajcc_stage', 'data_empty', 'IhdII NOS', 'Stage X', 'data_loss']},
                              'patient_gender': {'delete': ['data_loss', 'data_empty']},
                              'patient_race': 'ori',
                              'to_last_known_alive': 'mid',
                              't_stage': {'T1': ['T1b', 'T1a1', 'T1', 'T1b1', 'T1a', 'T1c', 'T1b2'],
                                            'T2': ['T2', 'T2a', 'T2a2', 'T2b', 'T2c', 'T2a1'],
                                            'T3': ['T3b', 'T3a', 'T3', 'T3c'],
                                            'T4': ['T4d', 'T4c', 'T4', 'T4b', 'T4a', 'T4e'],
                                            'delete': ['data_empty', 'data_loss', 'TX', 'T0', 'Tis']},
                              'n_stage': {'N0': ['N0', 'N0 (mol+)', 'N0 (i+)', 'N0 (i-)'],
                                          'N1': ['N1mi', 'N1a', 'N1', 'N1b', 'N1c'],
                                          'N2': ['N2a', 'N2', 'N2c', 'N2b'],
                                          'N3': ['N3a', 'N3b', 'N3', 'N3c'],
                                          'delete': ['data_empty', 'empty_data', 'data_loss', 'NX']},
                              'm_stage': {'M0': ['cM0 (i+)', 'M0'],
                                          'M1': ['M1a', 'M1', 'M1b', 'M1c'],
                                          'delete': ['MX', 'data_loss', 'data_empty']}
                              }
        self.clicfeat_dict['patient_age'] = {}
        self.clicfeat_dict['patient_weight'] = {}
        for i in range(20):
            self.clicfeat_dict['patient_age'][str(i * 5 + 1) + ' ~ ' + str(i * 5 + 5)] = \
                [str(j) for j in range(i * 5 + 1, i * 5 + 6)]
        self.clicfeat_dict['patient_age']['delete'] = [
            'data_loss', 'data_empty']
        for k in range(4, 12):
            self.clicfeat_dict['patient_weight'][str(k * 10 + 1) + ' ~ ' + str(k * 10 + 10)] = \
                [str(h) for h in range(k * 10 + 1, k * 10 + 11)]
        self.clicfeat_dict['patient_weight']['0 ~ 40'] = [
            str(l) for l in range(0, 41)]
        self.clicfeat_dict['patient_weight']['>120'] = [
            str(a) for a in range(120, 300)]
        self.clicfeat_dict['patient_weight']['delete'] = [
            'data_loss', 'data_empty']
        self.ext = 10
        self.istypical = True
        self.tpcl_mut = 'TP53'
        self.clicfeat_dict['ICD_O3_site'] = {'head_neck':    ['C00', 'C01', 'C02', 'C03', 'C04',
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
        self.classreplace = {'1111': 'Class_7', '1112': 'Class_8', '1120': 'Class_1', '1210': 'Class_2', '1220': 'Class_3',
                             '2110': 'Class_4', '2120': 'Class_5', '2211': 'Class_9', '2212': 'Class_10', '2220': 'Class_6'}
        self.class_list_costom = ['Class_1', 'Class_2', 'Class_3', 'Class_4',
                                  'Class_5', 'Class_6', 'Class_7', 'Class_8',
                                  'Class_9', 'Class_10']
