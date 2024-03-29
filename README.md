# Reconstruction-coreg-pipeline
Code for coregistration of intracranial electrode used by the Cortical Physiology Laboratory at Massachusetts General Hospital, PI: Sydney Cash, https://cashlab.mgh.harvard.edu/

This code uses a number of different open-source packages to localize electrodes implanted in the brain in the case of intracranial monitoring for the treatment of epilepsy. The idea is to use components of several packages that allows reconstruction and visualization of electrode locations to match to recorded intracranial neurophysiological data, whether in examining epileptiform activity, understanding cognitive processing, or examining the effects of direct electrical stimulation. 

This pipeline is focused on remaining in the native space of the patient as much as possible. For this reason, we do not use many components already provided by such software packages BrainStorm or LeadDBS. 

The design of the code is modular, in that certain products of the pipeline can be separated from others. Further, there are other options for the electrode localization (shown below).

## Useful code (either used in this pipeline or optional code):

MATLAB: https://www.mathworks.com/products/matlab.html

Freesurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall

iELVis: https://www.mathworks.com/matlabcentral/fileexchange/57317-ielvis

Fieldtrip: https://www.fieldtriptoolbox.org/ 

MMVT-lite: https://github.com/pelednoam/mmvt_lite

## For de-identification:
Fieldtrip: https://www.fieldtriptoolbox.org/ 

## For standardizing the data in iEEG BIDS format: 
BIDS starter kit WIKI: https://bids-standard.github.io/bids-starter-kit/

BIDS starter kit GitHub: https://github.com/bids-standard/bids-starter-kit 

## Optional software which would allow you to follow the same pipeline include the below:
Mango: https://ric.uthscsa.edu/mango/ (can be used in place of FreeView)

LeadDBS: https://www.lead-dbs.org/ 

LeGUI: https://github.com/Rolston-Lab/LeGUI


![Flowchart2](https://user-images.githubusercontent.com/11430978/204407584-68956d0a-a07f-426b-a9fa-fa0c4e88482d.png)

