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
data are extracted intermittently every 40 seconds (1 rest block, 1 task block) 
An similar approach has been used for
the localizer sessions whose data are used to generate the baseline stimuli
and all the files needed for the NF run.

This script requires a set of files already estimated using the utilities
of the rtRSA.



"""
from psychopy import core, visual, event
import os, pickle, serial, sys
import numpy as np
from expyriment_stash.extras.expyriment_io_extras import tbvnetworkinterface
from tkinter import filedialog
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
    
    line = ax.plot(rtRSAObj.RS_coords[0:1,0],rtRSAObj.RS_coords[0:1,1], color = 'dimgray')
    
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
    
    
    plt.savefig(os.path.join(outdir,'initial_img.png'), dpi=200,facecolor='dimgray', bbox_inches='tight')  
    
    # put it in an ImageStim
    image.setImage(os.path.join(outdir,'initial_img.png'))
    
    return scat, ax, line


def create_feedback(scat,ax, line, idx_fb,stimulus_positions,img,win, RScoords):
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
    
    elif idx_fb > 0 and idx_fb < 4:
        #from the second feedback to the fourth
        new_scat = scat
        #delete the reference to the point plotted for the previous iteration
        new_scat.set_offsets(np.delete(scat.get_offsets(), 0, axis=0))

        line.pop(0).remove()
        #plot the new position of the current mental state and, in dashed lines,
        #the trajectory until now        
        new_scat = ax.scatter(stimulus_positions[idx_fb,0],stimulus_positions[idx_fb,1], 
                    marker = '*',s=200, color = 'red', edgecolors='black',zorder=3)
        line = ax.plot(stimulus_positions[:idx_fb+1,0],stimulus_positions[:idx_fb+1,1], '--',
                        color = 'black',alpha=0.2,zorder=1) 

    else:
        #from the fourth feedback onwards
        new_scat = scat
        #delete the reference to the point plotted for the previous iteration
        new_scat.set_offsets(np.delete(scat.get_offsets(), 0, axis=0))
        
        line.pop(0).remove()
        #plot the new position of the current mental state and, in dashed lines,
        #the trajectory until now        
        new_scat = ax.scatter(stimulus_positions[idx_fb,0],stimulus_positions[idx_fb,1], 
                    marker = '*',s=200, color = 'red', edgecolors='black',zorder=3)
        line = ax.plot(stimulus_positions[(idx_fb-4):(idx_fb+1),0],stimulus_positions[(idx_fb-4):(idx_fb+1),1], '--',
                        color = 'black',alpha=0.2,zorder=1)

    ax.set_facecolor('dimgray')

	#scale plot
    
    xmin = np.min(np.append(RScoords[:,0], stimulus_positions[idx_fb,0]))
    xmax = np.max(np.append(RScoords[:,0], stimulus_positions[idx_fb,0]))
    ymin = np.min(np.append(RScoords[:,1], stimulus_positions[idx_fb,1]))
    ymax = np.max(np.append(RScoords[:,1], stimulus_positions[idx_fb,1]))
    dx = (xmax-xmin)*0.1
    dy = (ymax - ymin)*0.1
    ax.set_xlim(xmin - dx , xmax + dx)
    ax.set_ylim(ymin - dy, ymax + dy)
    

    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
            
    #save figure
    plt.savefig(os.path.join(outdir,'tvals_Trial' + str(idx_fb)+ '.png'),
               facecolor='dimgray', edgecolor='none', dpi=200, bbox_inches='tight')
    
    #load and display the saved figure
    image.setImage(os.path.join(outdir,'tvals_Trial' + str(idx_fb)+ '.png'))
    image.draw()
    win.flip()

    return new_scat, line





#%%###########################################################################
#                           Create an rtRSA object                           #
##############################################################################


#get current directory
wdir = os.getcwd()

#get the directory where there are the RSA data of the subject
sub_json_file = filedialog.askopenfilename(title='Select the config json file of the subject')

#empty RSA object
rtRSAObj = nfrsa.rtRSA(' ',2,' ')

#load properties in the just created rtRSAObj
#the config.json file can be created by using one of the class method
#the file si wirtten after the estimation of the RS
rtRSAObj.load(sub_json_file)

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
audio_stim_name = input('Insert the filename of the audio stimulus: '  )
nr_NF_session = input('Number of the NF session (1,2,3...): ')

#path of the stimulus
stimulus_path = os.path.join(wdir,'sounds/stimuli/'+audio_stim_name +'.wav')
stimulus = Sound(stimulus_path)
#path of the stop audio file
stop_wav= os.path.join(wdir,'sounds/stop.wav')
stop_stim = Sound(stop_wav)


#create a new output directory for the FB images
outdir = os.path.join(os.path.abspath(os.path.join(sub_json_file,os.pardir)),
                      audio_stim_name + '_output_' + str(nr_NF_session))
if not os.path.exists(outdir):
    os.makedirs(outdir)

#condition timings
baselines = pickle.load(open(os.path.join(wdir,'prt/long_intermittent/baselines.pkl'),'rb'))
tasks = pickle.load(open(os.path.join(wdir,'prt/long_intermittent/tasks.pkl'),'rb'))
feedbacks = pickle.load(open(os.path.join(wdir,'prt/long_intermittent/feedbacks.pkl'),'rb'))

nr_of_trials = feedbacks.shape[0]

#variable to store positions of the stimulus in time
stimulus_positions = np.empty((nr_of_trials+1,2))
feedback_distances = []
all_tvalues = []

#index of the feedbacks
idx_fb = 0

#index of the contrast map. TBV needs to know which contrast map we want
idx_ctr = 0

#%%############################################################################
#                             Reading data from TBV                           #
###############################################################################
timepoint_timing = []

scat, ax, line = create_baseline_figure(rtRSAObj,image)
print('Baseline figure created!')


#just display the first fixation
fixation.draw()
win.flip()

#just display the RS without feedback about the current brain pattern
# image.setImage(os.path.join(outdir,'initial_img.png')) 
# image.draw() 
# win.flip() 

real_run = 'n'
if  real_run == 'y':
    
    #These lines hold true only for a scanner with a trigger on a serial port!!
    #Please check your scan before use it!!!
    print('N.B. This works only for a scanner with a trigger based on serial port!!!')
    serFmri = serial.Serial('COM1', 57600)
    prevState = serFmri.getDSR()

    while serFmri.getDSR() == prevState:
        print('Waiting scanner....')

else:
#USEFUl FOR OTHER SCANNERS AND SIMULATIONS
    while '5' not in event.getKeys(['5']):
        print('Waiting scanner....')

#%% 


globalClock.reset()    
print("First trigger!")

#it waits until the first time point is processed by TBV to be sure to 
#read correct data from TBV settings file
CurrTimePoint = 0
while TBV.get_current_time_point()[0] < 1:
    NrOfTimePoints = TBV.get_expected_nr_of_time_points()[0]
    NrOfROIs = TBV.get_nr_of_rois()[0]
    #print('ROIs:',NrOfROIs)
    print('Waiting TBV....')

#boolean flag. It is True only if we have the coordinates of the ROI
tbv_coords_ok = False 


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
            
            if not tbv_coords_ok:
                
                #matching the current FMR coords with the ones of the localizer
                raw_nf_coords = np.array(TBV.get_all_coords_of_voxels_of_roi(0)[0])
                raw_nf_coords = raw_nf_coords[np.lexsort((raw_nf_coords[:,2], raw_nf_coords[:,1],raw_nf_coords[:,0]))]                    
                nf_coords = rtRSAObj.match_coords(np.array(raw_nf_coords))
                tbv_coords_ok = True #print(nf_coords)


            #Rest                        
            if CurrTimePoint in baselines[:,0]:
                fixation.draw()
                win.flip()

            
            #Imagination
            elif CurrTimePoint in tasks[:,0]:
                print('stimulus')
                stimulus.play()
                fixation.draw()
                win.flip()
                             
                
            #end of the task
            elif CurrTimePoint in tasks[:,1]:                            
                stop_stim.play()
                print('stop')
                
                
            elif CurrTimePoint == feedbacks[0,0]:
                image.setImage(os.path.join(outdir,'initial_img.png')) 
                image.draw() 
                win.flip() 

            #extract the map and plot the current position
            elif CurrTimePoint in feedbacks[1:,0]:
                
                #extracting tvalues from the ROI
                #in this experimental paradigm we have only one contrast
                tvalues = [TBV.get_map_value_of_voxel(idx_ctr,coords)[0] 
                                for coords in nf_coords]

                #estimate nwe stimulus coordinates
                stimulus_positions[idx_fb,0],stimulus_positions[idx_fb,1],tmp_dist = rtRSAObj.target_positioning(tvalues)
                
                #create the feedback
                if idx_fb == 0:
                    new_scat, line = create_feedback(scat,ax, line,idx_fb,stimulus_positions,image,win, rtRSAObj.RS_coords)
                else:
                    new_scat, line = create_feedback(new_scat,ax, line, idx_fb,stimulus_positions,image,win, rtRSAObj.RS_coords)
              
                #increment the index of the contrast map
                idx_fb += 1
                feedback_distances.append(tmp_dist)
                all_tvalues.append(tvalues)

            #pattern of the full time-series GLM
            elif CurrTimePoint == NrOfTimePoints:
                
                os.system('pause')
                #contrast number is fixed for the simulation
                #extract tvalues at the corresponding coordinates
                tvalues = [TBV.get_map_value_of_voxel(idx_ctr,coords)[0] 
                           for coords in nf_coords]
                #print(tvalues)
                #estimate new stimulus coordinates
                stimulus_positions[idx_fb,0],stimulus_positions[idx_fb,1],tmp_dist = rtRSAObj.target_positioning(tvalues)
                feedback_distances.append(tmp_dist)
                all_tvalues.append(tvalues)
                #create the feedback
                new_scat, line = create_feedback(new_scat,ax, line, idx_fb,stimulus_positions,image,win, rtRSAObj.RS_coords)
                
                print('Last time point!')

                with open(os.path.join(outdir,'timepoint_timing_incrementalGLM.txt'), 'w') as f:
                    for item in timepoint_timing:
                        f.write("%s\n" % item) 
                f.close()            
                
                pickle.dump(stimulus_positions,open(os.path.join(outdir,'feedback_positions.pkl'),'wb'))
                pickle.dump(feedback_distances,open(os.path.join(outdir,'feedback_distances.pkl'),'wb'))
                pickle.dump(all_tvalues,open(os.path.join(outdir,'all_tvalues.pkl'),'wb'))
                
                
                win.close()
                       
                print('Bye! See you soon!')
                sys.exit()
                    
              







    