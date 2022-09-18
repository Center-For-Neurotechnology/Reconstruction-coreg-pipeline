
clear 
addpath(genpath('X:\Projects\Lab_Materials\Analysis_Tools_and_Software\fieldtrip-20220202\'))

PatientName='sub-0d9l';
wd = ['Y:\StimDataBackup\Data_Stimulation\DeIdentifiedDataSet\'];

DirVal = fullfile(wd , ['derivatives\freesurfer\',PatientName], 'surf'); %the surf folder with the freesurfer output
RASExcelDirectory = fullfile(wd, PatientName, 'ses-postimp','ieeg' );
PreopMRIDirectory = fullfile(wd , [PatientName], 'ses-preimp', 'anat'); 
PreopMRIFile=dir([PreopMRIDirectory,'\*.nii']);

PostopMRIDirectory = fullfile(wd , [PatientName], 'ses-postimp', 'anat'); 
PostopMRIFile=dir([PostopMRIDirectory,'\*.nii']);


%%
clear ct mri mri_anon ct_anon

Desig=PatientName;

Desig
%%
mri = ft_read_mri([PreopMRIDirectory,PreopMRIFile(1).name]);
%
cfg = [];
mri_anon = ft_defacevolume(cfg, mri);
%
ft_write_mri([PreopMRIDirectory,'\',Desig,'_anon.nii'], mri_anon.anatomy, 'transform', mri_anon.transform, 'dataformat', 'nifti');

%%
ct = ft_read_mri([PostopMRIDirectory,PostopMRIFile(1).name]);
%
cfg = [];
ct_anon = ft_defacevolume(cfg, ct);
%
ft_write_mri([PostopMRIDirectory,'\',Desig,'_anonPostOPScan.nii'], ct_anon.anatomy, 'transform', ct_anon.transform, 'dataformat', 'nifti');
NUMCh
