from __future__ import absolute_import, division, print_function

from psychopy import core, visual,event
from psychopy.hardware.emulator import launchScan
import os
import numpy as np 
from psychopy.sound import Sound


stimulus_name = input('Name of the stimulus (cane, gatto, martello, sedia): \n')

wdir = os.getcwd()

##check this path


path_stimuli = os.path.join(wdir,'stimuli')

stop_wav= os.path.join(wdir,'stop.wav')

stimulus = Sound(os.path.join(path_stimuli,stimulus_name+'.wav'))
stop_stim = Sound(stop_wav)



win = visual.Window(fullscr = False, color='dimgray', screen = 1,
                    size=(1024,768), colorSpace = 'rgb255')


# settings for launchScan:
MR_settings = {
    'TR': 1,     # duration (sec) per whole-brain volume
    'volumes': 400,    # number of whole-brain 3D volumes per scanning run
    'sync': '5', # character to use as the sync timing event; assumed to come at start of a volume
    'skip': 0,       # number of volumes lacking a sync pulse at start of scan (for T1 stabilization)
    'sound': False    # in test mode: play a tone as a reminder of scanner noise
    }

duration = MR_settings['volumes'] * MR_settings['TR']#initialization of the duration of the experiment based on the MRI volumes

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

#creation of the duration variables
task_duration, baseline_duration = 20, 20

end_stimuli = 420 #20*11baseline + 20*10baseline
baselines_onset = np.arange(1,end_stimuli,40)
baselines_offset = baselines_onset + 19

task_onset = np.arange(21,end_stimuli-20,40)
task_offset = task_onset + 19


globalClock = core.Clock()#clock initialization


tot_vol = launchScan(win, MR_settings, globalClock=globalClock,esc_key='escape')

fixation.draw()
win.flip()

#while tot_vol < MR_settings['volumes']:
while globalClock.getTime()< duration:
    
    
    
    allKeys = event.getKeys()
    for key in allKeys:
        if key == MR_settings['sync']:
     
            #update the total volumes counter
            print('Current volume: ',str(tot_vol))
    
            #display the images according to the timing onset predefined
            if tot_vol in baselines_onset[1:]:
                print(globalClock.getTime())
                stop_stim.play()
                    
            if tot_vol in task_onset: 
                print(globalClock.getTime())
                stimulus.play()
    
        tot_vol += 1
            
      
        
win.close()
core.quit()