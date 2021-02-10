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

import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface
from matplotlib import pyplot as plt
import pickle 
import os



#%%############################################################################
#                               TBV  interface settings                       #
###############################################################################

#create an instance to access TBV via network plugin
TBV = tbvnetworkinterface.TbvNetworkInterface('localhost',55555)



nf_coords = TBV.get_all_coords_of_voxels_of_roi(0)[0]

n_ctr = TBV.get_nr_of_contrasts()[0]

tvalues_ml = np.empty((len(nf_coords),9))
#tvalues_agr = np.empty((len(nf_coords),9))

for ctr in range(9):
    tvalues_ml[:,ctr]= np.array([TBV.get_map_value_of_voxel(ctr,coords)[0] 
                for coords in nf_coords])

#%%

#pickle.dump(tvalues_agr,open(os.getcwd()+'/normal_prt.pkl','wb'))
pickle.dump(tvalues_ml,open(os.getcwd()+'/not_normal_prt.pkl','wb'))

#%%

tval_diff = (abs(tvalues_ml)-abs(tvalues_agr))/abs(tvalues_agr)#+np.array(tvalues_agr))
labels = ['pred1','pred2','pred3','pred4','pred5','pred6','pred6','pred7','pred8','pred9']

signs = np.sign(tvalues_agr)*np.sign(tvalues_ml)

plt.plot(tval_diff*signs*100)
plt.legend(labels,fontsize=20)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30) 
plt.tight_layout()


from scipy.stats import pearsonr

corr = []
for i in range(tval_diff.shape[1]): 
    corr.append(pearsonr(tvalues_agr[:,i],tvalues_ml[:,i]))