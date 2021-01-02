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
import pickle
import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface
import matplotlib
matplotlib.use('Agg') #to avoid display of the plot
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from psychopy.sound import Sound
from rtrsa import nfrsa


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
print('rtRSA ready!\n')
print('Nr of voxels: ',+ len(rtRSAObj.func_coords))
print('Base stimuli name:')
print(rtRSAObj.conditions)

#%%############################################################################
#                               TBV  interface settings                       #
###############################################################################

#create an instance to access TBV via network plugin
TBV = tbvnetworkinterface.TbvNetworkInterface('localhost',55555)

win = visual.Window(fullscr=False,color='gray',screen=0,
                    size=(1024,768),colorSpace='rgb255') 

#creation of the cue for the image
image = visual.ImageStim(win, 
                     image=None, 
                     mask=None, 
                     units='pix', 
                     pos=(0.0, 0.0), 
                     size=600, 
                     ori=0.0, 
                     color=(1.0, 1.0, 1.0), 
                     colorSpace='rgb', 
                     contrast=1.0, 
                     opacity=1.0, 
                     depth=0, 
                     interpolate=False, 
                     flipHoriz=False, 
                     flipVert=False, 
                     texRes=128, 
                     name=None, 
                     autoLog=None, 
                     maskParams=None)

#creation of the cue for the frame around the image (it is a text stimulus)
fixation = visual.TextStim(win, 
                           text='+', 
                           font='Arial', 
                           pos=(0.0, 0.0), 
                           depth=0, 
                           rgb=None, 
                           color='white', 
                           colorSpace='rgb', 
                           opacity=1.0, 
                           contrast=1.0, 
                           units='', 
                           ori=0.0, 
                           height=0.15, 
                           antialias=True, 
                           bold=False, 
                           italic=False, 
                           alignHoriz='center', 
                           alignVert='center', 
                           fontFiles=(), 
                           wrapWidth=None, 
                           flipHoriz=False, 
                           flipVert=False, 
                           languageStyle='LTR', 
                           name=None, 
                           autoLog=None)

#clock initialization
globalClock = core.Clock()

#%%############################################################################
#                      THESE VARIABLES DEPEND ON                             #
#                      THE EXPERIMENTAL DESIGN                               #
###############################################################################

#path of the stimulus
stimulus_path = os.path.join(wdir,'sounds/stimuli/cat.wav')
stimulus = Sound(stimulus_path)

stop_wav= os.path.join(wdir,'sounds/stop.wav')
stop_stim = Sound(stop_wav)

outdir = os.path.join(wdir,'output')

#condition timings
baselines = pickle.load(open(os.path.join(wdir,'prt/baselines.pkl'),'rb'))
tasks = pickle.load(open(os.path.join(wdir,'prt/tasks.pkl'),'rb'))
feedbacks = pickle.load(open(os.path.join(wdir,'prt/feedbacks.pkl'),'rb'))
fb_duration = feedbacks[0,1]-feedbacks[0,0] +1

nr_of_trials = feedbacks.shape[0]

#variable to store positions of the stimulus in time
stimulus_positions = np.empty((nr_of_trials+1,2))

#index of the constrast maps
idx_ctr = 0


#%%############################################################################
#                     Create a plot of the base stimuli                       #
###############################################################################

#setting a basic plot figure
plt.ion()
fig = plt.figure()

plt.figure(facecolor='dimgray')
ax = plt.axes()
ax.scatter(rtRSAObj.RS_coords[:,0],rtRSAObj.RS_coords[:,1],s=150, c='yellow',
           edgecolors = 'black')
# Setting the background color
ax.set_facecolor('dimgray')
plt.xticks([])
plt.yticks([])
plt.axis('off')
for label, x, y in zip(rtRSAObj.conditions, rtRSAObj.RS_coords[:,0],rtRSAObj.RS_coords[:,1]):
    plt.annotate(label, xy=(x, y), xytext=(20, -20),size=15,
                 textcoords='offset points', ha='right', va='bottom',
                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=1),
                 arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0')) 


yellow_circle = mlines.Line2D([], [], color='yellow',  linestyle='None', 
                              marker='o',markeredgecolor='black', markeredgewidth=0.5,
                          markersize=15, label='Base stimuli')
red_star = mlines.Line2D([], [], color='red',linestyle='None',  marker='*',
                         markeredgecolor='black', markeredgewidth=0.5,
                          markersize=15, label='Current stimulus')

plt.legend(handles=[yellow_circle, red_star],
    loc='upper center', bbox_to_anchor=(0.5, -0.05),ncol=4)




#%%############################################################################
#                             Reading data from TBV                           #
###############################################################################
timepoint_timing = []


#USEFUl FOR OTHER SCANNERS AND SIMULATIONS
while '5' not in event.getKeys(['5']):
    print('Waiting scanner....')


# serFmri = serial.Serial('COM1', 57600)
# prevState = serFmri.getDSR()

# while serFmri.getDSR() == prevState:
#     print('Waiting scanner....')
   

globalClock.reset()    
print("First trigger!")


#it waits until the first time point is processed by TBV to be sure to 
#read correct data from TBV settings file
CurrTimePoint = 0
while TBV.get_current_time_point()[0] < 1:
    NrOfTimePoints = TBV.get_expected_nr_of_time_points()[0]
    NrOfROIs = TBV.get_nr_of_rois()[0]
    print('Waiting TBV....')

raw_nf_coords = []

print("OK let's go! Expected TPs: " + str(NrOfTimePoints))    
#general loop
while TBV.get_current_time_point()[0] <= NrOfTimePoints:
    
    if CurrTimePoint !=  TBV.get_current_time_point()[0]:
        timepoint_timing.append([CurrTimePoint ,globalClock.getTime()])
        #update current timepoint
        CurrTimePoint = TBV.get_current_time_point()[0]
        print('Current time point:',str(CurrTimePoint))
        
            
        #looking for a ROI
        if NrOfROIs == 0:                
            print('Please add a ROI')
            NrOfROIs = TBV.get_nr_of_rois()[0]
            
            
        else:
            
            if not raw_nf_coords:
                raw_nf_coords = TBV.get_all_coords_of_voxels_of_roi(0)[0]
                    
                nf_coords = rtRSAObj.match_coords(raw_nf_coords)
            #needed to avoid accessing to timepoint -1 (fake) or timepoint 0
            if CurrTimePoint > 1 :
                #THE ACTUAL EXPERIMENT STARTS ONLY IF THERE IS A ROI!!!!!!#

                #coordinates of voxels of the VOI (index = 0)
                # raw_nf_coords = TBV.get_all_coords_of_voxels_of_roi(0)[0]
                
                # nf_coords = rtRSAObj.match_coords(raw_nf_coords)
                
                #setting the fixation for the baseline
                if CurrTimePoint in baselines[1:,0]:
                    fixation.draw()
                    win.flip()
                
                #showing only the frame for the imaginative task
                elif CurrTimePoint in tasks[:,0]:
                    print('stimulus')
                    stimulus.play()
                    fixation.draw()
                    win.flip()
                    
                elif CurrTimePoint in tasks[:,1]:
                    stop_stim.play()
                    print('stop')
                    
                #extract the map and plot the current position
                elif CurrTimePoint in feedbacks[:,0]:

                    #extractiing tvalues from the ROI
                    #in this experimental paradigm we have only one contrast
                    tvalues = [TBV.get_map_value_of_voxel(0,coords)[0] 
                                    for coords in nf_coords]

                    #estimate nwe stimulus coordinates
                    stimulus_positions[idx_ctr,:] = rtRSAObj.target_positioning(tvalues)
                    
                    
                    #plotting the new coordinates
                    if idx_ctr == 0:
                        print('Contrast:',idx_ctr+1)
                        plt.scatter(stimulus_positions[idx_ctr,0],
                                    stimulus_positions[idx_ctr,1], 
                                    marker = '*',s=200, color = 'red', 
                                    edgecolors='black')
                        ax.set_facecolor('dimgray')
                        plt.xticks([])
                        plt.yticks([])
                        plt.axis('off')
                    else:
                        print('Contrast:',idx_ctr)
                        plt.scatter(stimulus_positions[:idx_ctr,0],stimulus_positions[:idx_ctr,1], 
                                    marker = '*',s=200, color = 'darkgray')
                        plt.scatter(stimulus_positions[idx_ctr,0],stimulus_positions[idx_ctr,1], 
                                    marker = '*',s=200, color = 'red', edgecolors='black')
                        #plotting the trajectory
                        plt.plot(stimulus_positions[:idx_ctr+1,0],stimulus_positions[:idx_ctr+1,1], '-',
                                    color = 'green')
                        ax.set_facecolor('dimgray')
                        plt.xticks([])
                        plt.yticks([])
                        plt.axis('off')
                            
                    #save figure
                    plt.savefig(os.path.join(outdir,'tvals_Trial' + str(idx_ctr)+ '.png'),
                               facecolor='dimgray', edgecolor='none', dpi=300)
                    #show the figure()
                    image.setImage(os.path.join(outdir,'tvals_Trial' + str(idx_ctr)+ '.png'))
                    image.draw()
                    win.flip()
                    core.wait(fb_duration)
                    
                    #increment the index of the contrast map
                    idx_ctr += 1
                    
                elif CurrTimePoint == NrOfTimePoints:
                    print('Last time point!')
                    #contrast number is fixed for the simulation
                    #extract tvalues at the corresponding coordinates
                    tvalues = [TBV.get_map_value_of_voxel(0,coords)[0] 
                               for coords in nf_coords]
                    
                    #estimate new stimulus coordinates
                    stimulus_positions[idx_ctr,:] = rtRSAObj.target_positioning(tvalues)
                    
                    #plotting the new coordinates
                    print('Contrast:',idx_ctr+1)
                    plt.scatter(stimulus_positions[:idx_ctr,0],stimulus_positions[:idx_ctr,1], 
                                marker = '*',s=200, color = 'darkgray')
                    for label, x, y in zip(range(idx_ctr),stimulus_positions[:idx_ctr,0],stimulus_positions[:idx_ctr,1]):
                        plt.annotate(label,xy=(x, y), xytext=(-5, 5))
                    plt.scatter(stimulus_positions[idx_ctr,0],stimulus_positions[idx_ctr,1], 
                                marker = '*',s=200, color = 'red', edgecolors='black')
                    #plotting the trajectory
                    plt.plot(stimulus_positions[:idx_ctr+1,0],stimulus_positions[:idx_ctr+1,1],'-',
                                color = 'green')
                    ax.set_facecolor('dimgray')
                    plt.xticks([])
                    plt.yticks([])
                    plt.axis('off')

                            
                    #save figure
                    plt.savefig(os.path.join(outdir,'tvals_Trial' + str(idx_ctr)+ '.png'),
                               facecolor='dimgray', edgecolor='none', dpi=300)
                    #show the figure()
                    image.setImage(os.path.join(outdir,'tvals_Trial' + str(idx_ctr)+ '.png'))
                    image.draw()
                    win.flip()
                    core.wait(fb_duration)
                    
                    break

                else:
                    fixation.draw()
                    win.flip()
                       

                    
              
plt.close()

win.close()


with open(os.path.join(outdir,'timepoint_timing_incrementalGLM.txt'), 'w') as f:
    for item in timepoint_timing:
        f.write("%s\n" % item) 
f.close()            



    