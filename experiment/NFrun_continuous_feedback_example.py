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
import os,sys, pickle
import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface
#import matplotlib
#matplotlib.use('Agg') #to avoid display of the plot
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from psychopy.sound import Sound
from rtrsa import nfrsa

 

#%%############################################################################
#                                 FUNCTIONS                                   #
###############################################################################

    
def create_baseline_figure(rtRSAObj, image):
    '''
    
    Parameters
    ----------
    rtRSAObj : Object class rtRSA
        rtRSA object that contains all the info for the experiments (e.g. 
        coordinates of the RS space, the inversion matrix and the base stimuli).
    image : psychopy visual stimulus
        This object contains all the information to deliver a visual stimulus 
        to the subject.

    Returns
    -------
    scat : scatter plot instance
        Return the instance of the scatter plot created that has to be updated.
    ax : matplotlib axis
        Axis of the scatter plot created.

    '''
    
    fig,ax = plt.subplots()
    fig.set_facecolor('dimgray')
    scat = ax.scatter(rtRSAObj.RS_coords[:,0],rtRSAObj.RS_coords[:,1],s=150, c='yellow',
               edgecolors = 'black')
    ax.axis('off')
    
    for label, x, y in zip(rtRSAObj.conditions, rtRSAObj.RS_coords[:,0],rtRSAObj.RS_coords[:,1]):
        ax.annotate(label, xy=(x, y), xytext=(-20, 20),size=15,
                     textcoords='offset points', ha='right', va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=1),
                     arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'), zorder=1) 
    
    
    yellow_circle = mlines.Line2D([], [], color='yellow',  linestyle='None', 
                                  marker='o',markeredgecolor='black', markeredgewidth=0.5,
                              markersize=15, label='Base stimuli')
    red_star = mlines.Line2D([], [], color='red',linestyle='None',  marker='*',
                             markeredgecolor='black', markeredgewidth=0.5,
                              markersize=15, label='Current stimulus')
    
    ax.legend(handles=[yellow_circle, red_star],loc='upper center', bbox_to_anchor=(0.5, -0.05),ncol=4)
    
    
    plt.savefig(os.path.join(outdir,'initial_img.png'), dpi=200,facecolor='dimgray')  
    
    # put it in an ImageStim
    image.setImage(os.path.join(outdir,'initial_img.png'))
    
    return scat, ax


def create_feedback(scat,ax,idx_fb,stimulus_positions,img,win):
    '''
    

    Parameters
    ----------
    scat : scatter plot instance
        The initial scatter plot that has to be updated.
    ax : matplotlib axis
        The axes of the initial scatter plot.
    idx_fb : integer
        Index to acces the stimulus positions array.
    stimulus_positions : numpy array
        Array that contains the stimulus positions for each feedback display.
    curr_block : string
        It indicates if we are in a task or in a rest block.
    img : Psychopy visual object
        It is a PsychoPy instance that contain the image to deliver to the subject.
    win : PsychoPy windows
        It is the PsychoPy window that hosts all the visual stimuli.

    Returns
    -------
    new_scat : scatter plot instance
        The matplot instance of the new scatter plot. It will be used to delete
        in the next iteration the red star that indicates the previous feedback

    '''
      

    if idx_fb == 0:
        #first feedback        
        new_scat = ax.scatter(stimulus_positions[idx_fb,0],stimulus_positions[idx_fb,1], 
                    marker = '*',s=200, color = 'red', edgecolors='black')
        
    else:
        #from the second feedback onwards
        new_scat = scat
        #delete the reference to the point plotted for the previous iteration
        new_scat.set_offsets(np.delete(scat.get_offsets(), 0, axis=0))
        #plot the new position of the current mental state and, in dashed lines,
        #the trajectory until now        
        new_scat = ax.scatter(stimulus_positions[idx_fb,0],stimulus_positions[idx_fb,1], 
                    marker = '*',s=200, color = 'red', edgecolors='black',zorder=3)
        ax.plot(stimulus_positions[:idx_fb+1,0],stimulus_positions[:idx_fb+1,1], '--',
                        color = 'black',alpha=0.2,zorder=1)

    ax.set_facecolor('dimgray')
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
            
    #save figure
    plt.savefig(os.path.join(outdir,'tvals_Trial' + str(idx_fb)+ '.png'),
               facecolor='dimgray', edgecolor='none', dpi=200)
    
    #load and display the saved figure
    image.setImage(os.path.join(outdir,'tvals_Trial' + str(idx_fb)+ '.png'))
    image.draw()
    win.flip()

    return new_scat







#%%###########################################################################
#                           Create an rtRSA object                           #
##############################################################################


#get current directory
wdir = os.getcwd()

rtRSAObj = nfrsa.rtRSA('mician_itc',2,'euclidean')

#load properties in the just created rtRSAObj
#the config.json file can be created by using one of the class method
#the file si wirtten after the estimation of the RS
rtRSAObj.load(os.path.join(wdir,'mician_itc.json'))
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
audio_stim_name = input('Insert the filename of the audio stimulus:'  )

#path of the stimulus
stimulus_path = os.path.join(wdir,'sounds/stimuli/'+audio_stim_name +'.wav')
stimulus = Sound(stimulus_path)

stop_wav= os.path.join(wdir,'sounds/stop.wav')
stop_stim = Sound(stop_wav)

outdir = os.path.join(wdir,'output')

#condition timings
baselines = pickle.load(open(os.path.join(wdir,'prt/continuous/baselines.pkl'),'rb'))
tasks = pickle.load(open(os.path.join(wdir,'prt/continuous/tasks.pkl'),'rb'))
feedbacks = pickle.load(open(os.path.join(wdir,'prt/continuous/feedbacks.pkl'),'rb'))

nr_of_trials = feedbacks.shape[0]

#variable to store positions of the stimulus in time
stimulus_positions = np.empty((nr_of_trials+1,2))

#index of the feedbacks
idx_fb = 0

#index of the contrast map. TBV needs to know which contrast map we want
idx_ctr = 0

#%%############################################################################
#                             Reading data from TBV                           #
###############################################################################
timepoint_timing = []

scat, ax = create_baseline_figure(rtRSAObj,image)
print('Baseline figure created!')



#USEFUl FOR OTHER SCANNERS AND SIMULATIONS
while '5' not in event.getKeys(['5']):
    print('Waiting scanner....')


#serFmri = serial.Serial('COM1', 57600)
#prevState = serFmri.getDSR()

#while serFmri.getDSR() == prevState:
#    print('Waiting scanner....')
   

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
while TBV.get_current_time_point()[0] <= NrOfTimePoints+1:
    
    if CurrTimePoint !=  TBV.get_current_time_point()[0]:
        timepoint_timing.append([CurrTimePoint ,globalClock.getTime()])
        #update current timepoint
        CurrTimePoint = TBV.get_current_time_point()[0]
        print('Current time point:',str(CurrTimePoint))
        
            
        #looking for a ROI
        if NrOfROIs == 0:                
            print('Please add a ROI')
            NrOfROIs = TBV.get_nr_of_rois()[0]
            #THE ACTUAL EXPERIMENT STARTS ONLY IF THERE IS A ROI!!!!!!#
                       
        else:
            
            if not raw_nf_coords:
                
                #matching the current FMR coords with the ones of the localizer
                raw_nf_coords = TBV.get_all_coords_of_voxels_of_roi(0)[0]                    
                nf_coords = rtRSAObj.match_coords(np.array(raw_nf_coords))
                
            #needed to avoid accessing to timepoint -1 (fake) or timepoint 0
            # while CurrTimePoint < baselines[0,0] :
            #     fixation.draw()
            #     win.flip()
               
            #showing only the frame for the imaginative task
            if CurrTimePoint in tasks[:,0]:
                print('stimulus')
                stimulus.play()
                                    
            elif CurrTimePoint in baselines[:,0]:
                stop_stim.play()
                print('stop')
                fixation.draw()
                win.flip()
                
            #extract the map and plot the current position
            elif CurrTimePoint in feedbacks:
                
                #extracting tvalues from the ROI
                #in this experimental paradigm we have only one contrast
                tvalues = [TBV.get_map_value_of_voxel(idx_ctr,coords)[0] 
                                for coords in nf_coords]

                #estimate nwe stimulus coordinates
                stimulus_positions[idx_fb,:] = rtRSAObj.target_positioning(tvalues)
                
                #create the feedback
                if idx_fb == 0:
                    new_scat = create_feedback(scat,ax,idx_fb,stimulus_positions,image,win)
                else:
                    new_scat = create_feedback(new_scat,ax,idx_fb,stimulus_positions,image,win)
              
                #increment the index of the contrast map
                idx_fb += 1
                
            elif CurrTimePoint == NrOfTimePoints:
                #contrast number is fixed for the simulation
                #extract tvalues at the corresponding coordinates
                tvalues = [TBV.get_map_value_of_voxel(idx_ctr,coords)[0] 
                           for coords in nf_coords]
                print(tvalues)
                #estimate new stimulus coordinates
                stimulus_positions[idx_fb,:] = rtRSAObj.target_positioning(tvalues)
                
                #create the feedback
                new_scat = create_feedback(new_scat,ax,idx_fb,stimulus_positions,image,win)
                
                print('Last time point!')

                break

                       

                    
              
#plt.close()



with open(os.path.join(outdir,'timepoint_timing_incrementalGLM.txt'), 'w') as f:
    for item in timepoint_timing:
        f.write("%s\n" % item) 
f.close()            

pickle.dump(stimulus_positions,open(os.path.join(outdir,'feedback_positions.pkl'),'wb'))


win.close()



    