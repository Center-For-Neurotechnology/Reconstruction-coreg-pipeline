%% Sets up the directories and the files to de-identify/deface
clear 
addpath(genpath('X:\Projects\Lab_Materials\Analysis_Tools_and_Software\fieldtrip-20220202\'))

PatientName='sub-5o1r';
wd = ['Y:\ReconPipelinePaper\Data\'];

DirVal = fullfile(wd , ['derivatives\freesurfer\',PatientName], 'surf'); %the surf folder with the freesurfer output
RASExcelDirectory = fullfile(wd, PatientName, 'ses-postimp','ieeg' );
PreopMRIDirectory = fullfile(wd , [PatientName], 'ses-preimp', 'anat'); 
PreopMRIFile=dir([PreopMRIDirectory,'\*.nii']);

PostopMRIDirectory = fullfile(wd , [PatientName], 'ses-postimp', 'anat'); 
PostopMRIFile=dir([PostopMRIDirectory,'\*.nii']);

clear post_anon mri mri_anon post

Desig=PatientName;
Desig
%% Reads the preoperative MRI
mri = ft_read_mri([PreopMRIDirectory,'\',PreopMRIFile(1).name]);
%
cfg = [];
mri_anon = ft_defacevolume(cfg, mri);
%
ft_write_mri([PreopMRIDirectory,'\',Desig,'_anon.nii'], mri_anon.anatomy, 'transform', mri_anon.transform, 'dataformat', 'nifti');

%% Reads the postoperative scan (CT or MRI)
post = ft_read_mri([PostopMRIDirectory,'\',PostopMRIFile(1).name]);
%
cfg = [];
post_anon = ft_defacevolume(cfg, post);
%
ft_write_mri([PostopMRIDirectory,'\',Desig,'_anonPostOPScan.nii'], post_anon.anatomy, 'transform', post_anon.transform, 'dataformat', 'nifti');
