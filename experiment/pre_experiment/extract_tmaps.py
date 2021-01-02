# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:20:30 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

"""


from rtrsa import utils
from tkinter import filedialog



ip = input('Insert TBV IP (default: localhost):  ')

if not ip:
    ip='localhost'

voi = int(input('Insert VOI index (default is 0):   '))

ctr = int(input('Insert contrast index (default is 0):   '))

single_maps = filedialog.askdirectory(title='Select a directory for the storing of the t-values maps')

out_name = input('Insert name of the ouput:   ')

utils.TBV_value_extractor(ip,voi,ctr,single_maps,out_name)