# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:24:48 2019

@author: Admin
"""
#from __main__ import *

import os
import importlib
import time
import Binning_optimization_v3 
import subprocess


DegFileimport = 'Strat_HPRD'
min_sep = 3
MDgenes ='90pc_S1.txt'

DegFname = DegFileimport + '_degree.txt'
ans0 = Binning_optimization_v3.calculateDegree('HPRD_GC_PPI',DegFileimport)

ans=Binning_optimization_v3.Binning(DegFname, min_sep)
ans2 = Binning_optimization_v3.MDStratification(DegFname, MDgenes)


    
