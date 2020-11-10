# ====================================
# Used for get stat data
# Written by Ji Hongchen
# 20200216
# ====================================

import matplotlib.pyplot as plt
import SATA_PRETREAT
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
import pandas as pd
import numpy as np
import scipy
import matplotlib as mlp
mlp.use('MacOSX')


class SataMethod():
    """
    paramater:
        |-> dataframe:      the dataframe output by test.py
        |-> class_col:      the class column pretreaded by SATA_PRETREAT.py
    """

    def __init__(self, dataframe,
                 class_col='class_result'):
        self.dataframe = dataframe
        self.class_col = class_col

    def _DropEmptyCol(self, drop_col=[]):
        """
        paramater:
            |-> drop_col: a list of columns that needed to delete empty
        """
        new_frame = self.dataframe
        for c in drop_col:
            new_frame = SATA_PRETREAT.DropEmptyData(new_frame,
                                                    column=c)
        return new_frame

    def _Neumeric(self, neu_col='', rep_zero=False):
        _neu_frame = self.dataframe
        _neu_frame[neu_col] = pd.to_numeric(
            _neu_frame[neu_col], errors='coerce')
        if rep_zero == True:
            _neu_frame[neu_col][_neu_frame[neu_col] < 0] = 0
        _neu_frame.fillna(0)
        return _neu_frame

    def Describe(self, descb_col=''):
        classes = np.unique(np.array(self.dataframe[self.class_col]))
        #_sata_frame = self._DropEmptyCol(drop_col=[descb_col])
        _sata_frame = self._Neumeric(neu_col=descb_col)
        #_sata_frame = self._Neumeric(neu_col=descb_col)
        index_dict = ['count', 'min', 'max', 'q25%',
                      'q50%', 'q75%', 'mean', 'var', 'std']
        describe_frame = pd.DataFrame(columns=classes, index=index_dict)
        for c in classes:
            _test_column = _sata_frame[_sata_frame[self.class_col]
                                       == c][descb_col]
            describe_frame[c]['count'] = _test_column.count()
            describe_frame[c]['max'] = _test_column.max()
            describe_frame[c]['min'] = _test_column.min()
            describe_frame[c]['q25%'] = _test_column.quantile(0.25)
            describe_frame[c]['q50%'] = _test_column.median()
            describe_frame[c]['q75%'] = _test_column.quantile(0.75)
            describe_frame[c]['mean'] = _test_column.mean()
            describe_frame[c]['var'] = _test_column.var()
            describe_frame[c]['std'] = _test_column.std()
        return describe_frame

    def CalTTest(self, describeframe):
        #_sata_frame = self._DropEmptyCol(drop_col=[ttest_col])
        t_test_frame = pd.DataFrame(columns=['T', 'P'])
        classes = describeframe.columns
        for i in range(len(classes)):
            for j in range(i+1, len(classes)):
                _mean_1 = describeframe[classes[i]]['mean']
                _mean_2 = describeframe[classes[j]]['mean']
                _std_1 = describeframe[classes[i]]['std']
                _std_2 = describeframe[classes[j]]['std']
                _n_1 = describeframe[classes[i]]['count']
                _n_2 = describeframe[classes[j]]['count']
                ttest_result = scipy.stats.ttest_ind_from_stats(mean1=_mean_1,
                                                                std1=_std_1,
                                                                nobs1=_n_1,
                                                                mean2=_mean_2,
                                                                std2=_std_2,
                                                                nobs2=_n_2)
                t_test_frame.loc[str(classes[i])+'_vs_' +
                                 str(classes[j])] = ttest_result
        return t_test_frame

    def ChiSeqTest(self, chisq_col='', recomb=[]):
        _sata_frame = self._DropEmptyCol(drop_col=[chisq_col])
        chiseq_frame = SATA_PRETREAT.ChiSquarePretreat(self.dataframe,
                                                       column=chisq_col,
                                                       class_column=self.class_col)
        chiseq_frame = chiseq_frame.fillna(0)
        if recomb == []:
            chiseq_array = np.array(chiseq_frame)
        else:
            chiseq_array = np.array(chiseq_frame[recomb])
        print(chiseq_array)
        chiseq_value = scipy.stats.chi2_contingency(chiseq_array)
        chiseq_frame['chi2'] = chiseq_value[0]
        chiseq_frame['P'] = chiseq_value[1]
        chiseq_frame['df'] = chiseq_value[2]
        return chiseq_frame

    def LogRankTest(self, outpath=''):
        #_sata_frame = self._DropEmptyCol(drop_col=['to_last_known_alive'])
        _sata_frame = self.dataframe
        classes = np.unique(np.array(self.dataframe[self.class_col]))
        km = KaplanMeierFitter()
        #fig = plt.figure()
        lr_frame = pd.DataFrame(columns=['stat_value', 'lr_p'])
        for i in range(len(classes)):
            for j in range(i+1, len(classes)):
                _T_1 = _sata_frame[_sata_frame[self.class_col]
                                   == classes[i]]['to_last_known_alive']
                _T_1 = pd.to_numeric(_T_1, errors='coerce')
                _E_1 = _sata_frame[_sata_frame[self.class_col]
                                   == classes[i]]['patient_status']
                _E_1 = pd.to_numeric(_E_1, errors='coerce')
                _T_2 = _sata_frame[_sata_frame[self.class_col]
                                   == classes[j]]['to_last_known_alive']
                _T_2 = pd.to_numeric(_T_2, errors='coerce')
                _E_2 = _sata_frame[_sata_frame[self.class_col]
                                   == classes[j]]['patient_status']
                _E_2 = pd.to_numeric(_E_2, errors='coerce')
                lr_p = logrank_test(_T_1, _T_2, _E_1, _E_2, alpha=0.95)
                lr_frame.loc[str(classes[i])+'__vs__'+str(classes[j])
                             ] = [lr_p.test_statistic, lr_p.p_value]
            try:
                km.fit(_T_1, event_observed=_E_1,
                       label='Class'+str(classes[i]))
            except UnboundLocalError:
                pass
            # km.plot()
            # km.fit(_T_2, event_observed=_E_2,
            #       label='Class'+str(classes[j]))
            try:
                km.plot()
            except AttributeError:
                pass
        if outpath != '':
            plt.savefig(outpath + '.png')
            plt.close()
            # plt.show()
        return (lr_frame)
