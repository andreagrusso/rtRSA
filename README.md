# rtRSA
A Python implementation of a real-time Representational Similarity Analysis for fMRI Neurofeedback experiment using Turbo-BrainVoyager.

If you use this tool please cite: (https://www.biorxiv.org/content/10.1101/2020.11.09.374397v2 in update...)

# How to use rt-SA 

  1) For a set of N stimuli/runs extract from Turbo-BrainVoyager the t-statistics relative to a GLM contrast from a single ROI by using extract_tmaps.py
  2) Create a rt-RSA object by using the extracted t-statistics. A rt-RSA object is defined by its *.json* file stored in the corresponding folder
  3) To run an experiment with the use of a rt-RSA object you need to load the *.json* file 

# Examples
Examples of possible Python scripts to run a rt-fMRI-NF experiment with the rt-RSA and different paradigm (i.e. continuous and intermeittent) are the following files:
  

  1) NFrun_int_FB_example.py  (feedback display after the task block)
  2) NFrun_long_int_FB_example.py (feedback display after the baseline block)
  3) NFrun_7T_paradigm_example.py (similar to nr. 2)
  4) NFrun_cnt_FB_example.py (feedback updates every 2s during the task block)
 
N.B. All examples are based on Turbo-BrainVoyager and its network plugin (https://www.brainvoyager.com/downloads/install_turbobrainvoyager.html) 

# Python packages needed
numpy_indexed; 
expyriment;
expyriment-stash;
numpy;
scipy; 
scikit-learn;
PsychoPy3 (to use the example experiments)



