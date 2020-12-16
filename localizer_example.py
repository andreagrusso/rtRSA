# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:56:46 2019

@author: Andrea Gerardo Russo
Biomedical Engineer, PhD candidate in Neuroscience
University of Salerno, Italy
e-mail: andrusso@unisa.it


This is a script to deliver visual stimuli (images in natural contest) to 
a subject during a fMRI acquisition.  

The 4 stimuli will be showed at least five times to increase the activation and
a block design paradigm will be used. Each block is composed of three parts:
    1) 2 seconds of visual presentation of the stimulus. Each image will be 
        presented in a squared white frame
    2) 20 seconds in which only the white frame will be displayed. Subject 
        would be asked to mentally imagine the stimulus showed before 
    3) 20 seconds of rest. This will serve as control condition to extract the
        t-maps
        


"""

from psychopy import core, visual,event
#from psychopy.hardware.emulator import launchScan
import os
import numpy as np 
from psychopy.sound import Sound
import serial



###############################################################################
#                              MRI SETTING                                    #          
###############################################################################

#######################Only For Skyra##########################################
#Getting Siemens Skyra trigger (it is based on voltage change)
#serFmri = serial.Serial('COM1', 57600)
#prevState = False
###############################################################################

stimulus_name = input('Name of the stimulus (cat, chair, dog, hammer): \n')




win = visual.Window(fullscr = True, color='dimgray', screen = 1,
                    size=(1024,768), colorSpace = 'rgb255')
print('Window size: ', win.size)



wdir = os.getcwd()

##check this path
outdir = os.path.join(wdir,"output") 

if not os.path.exists(outdir):
    os.mkdir(outdir)

path_stimuli = os.path.join(wdir,'sounds/stimuli')
stop_wav= os.path.join(wdir,'sounds/stop.wav')
    
if not os.path.exists(path_stimuli):
    print('Your sound stimuli are missing!\n Please add your .wav files to continue')
    
    exit()


###############################################################################
#                               TESTING SETTING                              #
###############################################################################
# settings for launchScan:
MR_settings = {
    'TR': 1,     # duration (sec) per whole-brain volume
    'volumes': 429,    # number of whole-brain 3D volumes per scanning run
    'sync': '5', # character to use as the sync timing event; assumed to come at start of a volume
    'skip': 0,       # number of volumes lacking a sync pulse at start of scan (for T1 stabilization)
    'sound': False    # in test mode: play a tone as a reminder of scanner noise
    }

duration = MR_settings['volumes'] * MR_settings['TR']#initialization of the duration of the experiment based on the MRI volumes

#Skyra in Salerno has a trigger that is based on the serial port! 
#We cannot use the launchScan()
serFmri = serial.Serial('COM1', 57600)
prevState = False



globalClock = core.Clock()#clock initialization
###############################################################################
#                               STIMULUS SETTINGS                             #
###############################################################################

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




#instantiating the sound stimulus to save time
stimulus = Sound(os.path.join(path_stimuli,stimulus_name+'.wav'))
stop_stim = Sound(stop_wav)
#initializing the index to get the first image in the array
stim_idx = 0

#creation of the duration variables
task_duration, baseline_duration = 20, 20

end_stimuli = 420 #20*11baseline + 20*10baseline
baselines_onset = np.arange(1,end_stimuli,40)
baselines_offset = baselines_onset + 19

task_onset = np.arange(21,end_stimuli-20,40)
task_offset = task_onset + 19

###############################################################################
#                                 USEFUL OUTPUTS                              #
###############################################################################

#timing of the stimuli
audio_timing = []


###############################################################################
#                               EMULATION                                     #
###############################################################################

# note: globalClock has been reset to 0.0 by launchScan()
#message = 'As you hear the word ' + stimulus_name.upper() + ' please start imaging the corresponding \
#image as clearly as possibile until you hear the word STOP. \
#Please keep your eyes on the fixation cross during the whole session.' 
#tot_vol = launchScan(win, MR_settings, globalClock=globalClock, wait_msg=message)

tot_vol = 0 #inital volume 
prevState = serFmri.getDSR()
           
while tot_vol<MR_settings['volumes']:
    
    currentState = serFmri.getDSR()
    if(currentState != prevState):
        if tot_vol == 0:
            globalClock.reset() #reset time at first pulse
        prevState = currentState
        tot_vol = tot_vol + 1
        
        
        win.flip()
        fixation.draw()
        #update the total volumes counter
        print('Current volume: ',str(tot_vol))
    
        #display the images according to the timing onset predefined
        if tot_vol in baselines_onset[1:]:
            print(globalClock.getTime())
            stop_stim.play()
                
        if tot_vol in task_onset: 
            audio_timing.append(globalClock.getTime())
            stimulus.play()
        
        
win.close()
 


with open(os.path.join(outdir,'audio_stim_timing_NF1.txt'), 'w') as f:
    for item in audio_timing:
        f.write("%s\n" % item) #write the complete item list (with brackets)
f.close()