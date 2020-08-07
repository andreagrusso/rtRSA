# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:54:25 2020

@author: Andrea Gerardo Russo, Biomedical Engineer
PhD candidate in Neuroscience
University of Salerno, Fisciano, Italy

This script provides an example on how to use the rtRSA approach in a 
rt-fMRI-NF experiment. As the rtRSA is built to work with Turbo-BrainVoyager
you need to have it installed on ypur laptop.

The experimental framework is built using PsychoPy3, although the rtRSA can
be used with any other stimulation packages. It is a block design with an 
imagery task. An auditory cue is delivered to the subject that is requested
to imagine the corresponding object until a stop cue is delivered. The brain 
data are extracted every 60 seconds (2 baseline,1 task) and provided as visual 
feedback to the subject for 5 seconds. An identical approach has been used for
the localizer sessions whose data are used to generate the baseline stimuli
and all the files needed for the NF run.

This script requires a set of files already estimated using the utilities
of the rtRSA.



"""
from psychopy import core, visual, event
import os
import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface
import matplotlib
matplotlib.use('Agg') #to avoid display of the plot
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from psychopy.sound import Sound
import json
from rtrsa import nfrsa, utils


#%%

#get current directory
wdir = os.getcwd()



#%%###########################################################################
#                           Create an rtRSA object                           #
##############################################################################

rtRSAObj = nfrsa.rtRSA('test',2,'pearson')

#load properties in the just created rtRSAObj
#the config.json file can be created by using one of the class method
#the file si wirtten after the estimation of the RS
rtRSAObj.load(os.path.join(wdir,'config.json'))


#%%############################################################################
#                               TBV  interface settings                       #
###############################################################################

#create an instance to access TBV via network plugin
TBV = tbvnetworkinterface.TbvNetworkInterface('localhost',55555)

    