%% Template Matlab script to create an BIDS compatible _electrodes.json file
% For BIDS-iEEG
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% Writing json files relies on the JSONio library
% https://github.com/gllmflndn/JSONio
% Make sure it is in the matab/octave path
%
% DHermes, 2017
% modified Jaap van der Aar 30.11.18

%%
% clear;
% % root_dir = ['Y:\StimDataBackup\Data_Stimulation\MG141\NetworkWake1\FT_Analyzed\UpdatedFT_Analyzed\bidsiEEG2\'];
% ieeg_project = 'templates';
% ieeg_sub = '01';
% ieeg_ses = '01';
% MRIDataAlignedTo='test';


electrodes_json_name = fullfile(root_dir, project_label, ...
                                ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
                                ['sub-' ieeg_sub ...
                                 '_ses-' ieeg_ses ...
                                 '_coordsystem.json']);

%% Required fields

loc_json.iEEGCoordinateSystem = 'ACPC'; % Defines the coordinate system for the iEEG electrodes.
% For example, "ACPC". See Appendix VIII: preferred names of Coordinate systems.
% If "Other" (for example: individual subject MRI), provide definition of the coordinate system
% in iEEGCoordinateSystemDescription.
% If positions correspond to pixel indices in a 2D image (of either a volume-rendering,
% surface-rendering, operative photo, or operative drawing), this must be "pixels".
% See section 3.4.1: Electrode locations for more information on electrode locations.

loc_json.iEEGCoordinateUnits = 'mm'; % Units of the _electrodes.tsv, MUST be "m", "mm", "cm" or "pixels".

%% Recommended fields

loc_json.iEEGCoordinateProcessingDescripton = 'Coordinate system with the origin at anterior commissure (AC), negative y-axis going through the posterior commissure (PC), z-axis going to a mid-hemisperic point which lies superior to the AC-PC line, x-axis going to the right'; 
% Freeform text description or link to document
% describing the iEEG coordinate system system in detail (for example: "Coordinate system with the origin
% at anterior commissure (AC), negative y-axis going through the posterior commissure (PC), z-axis
% going to a mid-hemisperic point which lies superior to the AC-PC line, x-axis going to the right")

loc_json.IntendedFor = [postimMRIName]; % This can be an MRI/CT or a file containing the operative photo, x-ray
% or drawing with path relative to the project folder. If only a surface reconstruction is available,
% this should point to the surface reconstruction file. Note that this file should have the same coordinate
% system specified in iEEGCoordinateSystem. (for example: "sub-<label>/ses-<label>/anat/sub-01_T1w.nii.gz")
% for example
% T1: "/sub-<label>/ses-<label>/anat/sub-01_T1w.nii.gz"
% Surface: "/derivatives/surfaces/sub-<label>/ses-<label>/anat/sub-01_T1w_pial.R.surf.gii"
% Operative photo: "/sub-<label>/ses-<label>/ieeg/sub-0001_ses-01_acq-photo1_photo.jpg"
% Talairach: "/derivatives/surfaces/sub-Talairach/ses-01/anat/sub-Talairach_T1w_pial.R.surf.gii"

loc_json.iEEGCoordinateProcessingDescription = ['Following implant, the preoperative T1-weighted MRI was aligned with a postoperative CT or MRI using volumetric image coregistration procedures and FreeSurfer scripts (http://surfer.nmr.mgh.harvard.edu). ',...
    'Electrode coordinates were manually determined from the CT or postoperative MRI in the patients’ native space and mapped using an electrode labeling algorithm (ELA) that registered each contact to a standardized cortical map. ',...
    'Surface projection for grid and strip electrodes was done using the snapgrid code in iElvis (https://github.com/iELVis/iELVis)']; % Has any projection been done on the electrode positions
% (for example: "surface_projection", "none").

loc_json.iEEGCoordinateProcessingReference = ['Original reconstructions done using Dykstra AR, Chan AM, Quinn BT, Zepeda R, Keller CJ, Cormier J, Madsen JR, Eskandar EN, Cash SS (2012) Individualized localization and cortical surface-based registration of intracranial electrodes. Neuroimage 59:3563–3570. ',...
'code in iElvis (https://github.com/iELVis/iELVis) was used for projecting the grid or strip electrodes ',...
'Further mapping was done using MMVT (https://github.com/pelednoam/mmvt) and mapped using an electrode labeling algorithm (ELA; https://github.com/pelednoam/ieil)']; % A reference to a paper that defines in more detail'}
% the method used to project or localize the electrodes

%% Write
jsonSaveDir = fileparts(electrodes_json_name);
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist, create: %s \n', jsonSaveDir);
end

json_options.indent = ' '; % this just makes the json file look prettier
% when opened in a text editor

try
    jsonwrite(electrodes_json_name, loc_json, json_options);
catch
    warning('%s\n%s\n%s\n%s', ...
            'Writing the JSON file seems to have failed.', ...
            'Make sure that the following library is in the matlab/octave path:', ...
            'https://github.com/gllmflndn/JSONio');
end