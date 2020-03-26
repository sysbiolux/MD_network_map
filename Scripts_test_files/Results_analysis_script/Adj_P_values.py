# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 14:05:43 2018

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 15:09:30 2018
Compute adj-p values 
Input: Csv file obtained from Analysis.py

@author: Admin
"""

import statsmodels.stats.multitest as smm
import pandas as pd



readcsv = pd.read_csv('Pvals_TEST_1_03.csv', header = None, names = ['Gene','No_of_Hits','Total_Hits','P_value','Av_cent','MD_cent','Rel_cent'])
Selrange = readcsv.loc[readcsv['Av_cent'] != 0]
Pvals_col = Selrange.iloc[:,3]      
pvalsForqvals = Pvals_col.apply(pd.to_numeric)     
[rej,pval_corr] = smm.fdrcorrection(pvalsForqvals,alpha=0.05)            
[rejMD, pval_corrMD,a,b] = smm.multipletests(pvalsForqvals, alpha=0.05, method='fdr_bh',is_sorted = False)
newcol = 'Q_value'
Qvalues = pd.DataFrame(pval_corr)
Qvalues.columns = ['Qvalues']
Selrange = Selrange.reset_index(drop=True)
Qvalues = Qvalues.reset_index(drop=True)
Alldata = Selrange.join(Qvalues)
Alldata.to_csv('P_Q_vals_HPRD_5000', sep='\t')            
