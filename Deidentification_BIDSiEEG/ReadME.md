# Cash lab iEEG BIDS formatting

Code in in this folder is to de-identify scans (by manual defacing) and save channel information and scans into the correct iEEG BIDS formatting for shared data. 

The example code where most the scripts are run together is "createBIDS_ieeg_jsoniEEGChannelsAll2.m". 

The code is modified from the bids-starter-kit:
https://github.com/bids-standard/bids-starter-kit/tree/main/matlabCode/ieeg

Following instructions for the iEEG BIDS format:

https://www.nature.com/articles/s41597-019-0105-7 

and here:

Brain Imaging Data Structure 1.7.0 Intacranial electroencephalography:

https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/04-intracranial-electroencephalography.html 

iEEG BIDS tutorial: 

https://bids-standard.github.io/bids-starter-kit/tutorials/ieeg.html 

Dependencies for this MATLAB pipeline include:

JSONio: a MATLAB/Octave JSON library

https://github.com/gllmflndn/JSONio 

Fieldtrip toolbox: 

https://www.fieldtriptoolbox.org/

https://github.com/fieldtrip/fieldtrip

Some toolboxes that help:

bids-MATLAB:

https://github.com/bids-standard/bids-matlab


