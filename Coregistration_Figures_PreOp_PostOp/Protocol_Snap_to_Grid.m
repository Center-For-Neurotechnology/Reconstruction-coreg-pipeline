%% Modular Reconstruction and Co-registration of Imaging from Implanted ECoG and SEEG Electrodes
%
% Code for Snapping Grids and Strips to Pial surface
%
% This version was completed on November 29th, 2022
% 
% For questions, please post issues to the GitHub
% 

Pt = 'sub-0t3i';
reconPath = '/MyPC/Documents/Recons'; 
[num,txt,all]=xlsread(fullfile(reconPath, Pt, [Pt '_RAS.xlsx'] ) );

%Number of RAS Coords
RAS_coords=[];
RAS_coords=num(1:length(num),:);

RAS_labels={};
RAS_labels=txt(1:length(txt),1);

% Adds the customized snap function from GitHub and the iELVis folder
addpath(genpath('/MyPC/Documents/MATLAB/Recons/Coregistration_Figures_PreOp_PostOp'))
addpath(genpath('/MyPC/Documents/MATLAB/Recons/CoregCode/iELVis'))

%% Separate each of the grids and strips into different variables
%G = RAS_coords(1:56,:);%Lateral Temporal Grid
IHA = RAS_coords(1:6,:);
IHP = RAS_coords(15:18,:);

GS = RAS_coords(19:34,:);
GI1 = RAS_coords(35:42,:);
GI2 = RAS_coords(43:50,:);

%Depths -- Will have to be adjusted manually below if you wish to move them
depths = RAS_coords(7:14,:);


%% Code for adjusting the locations of the grids and strips
surface_path = fullfile( reconPath , Pt, [Pt '_SurferOutput'], 'surf' );
side = 'r';
% surftype = 'pial-outer-smoothed';
[surf.vert surf.tri]=read_surf( fullfile( surface_path, [side  'h.pial-outer-smoothed'] ) );

% We recommend running each of the following lines one at a time
[gr1_coor_snapped] = snap2dural_energy_customized(IHA,surf);
[gr2_coor_snapped] = snap2dural_energy_customized(IHP,surf);
[gr3_coor_snapped] = snap2dural_energy_customized(GS,surf);
[gr4_coor_snapped] = snap2dural_energy_customized(GI1,surf);
[gr5_coor_snapped] = snap2dural_energy_customized(GI2,surf);

% Example advanced Usage
%[REF_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(REF,surf,nchoosek(1:length(REF),2));

%% Compile new coordinates
RAS_coords_new = [gr1_coor_snapped; gr2_coor_snapped; depths; gr3_coor_snapped; gr4_coor_snapped; gr5_coor_snapped];

R = RAS_coords_new(:,1);
A = RAS_coords_new(:,2);
S = RAS_coords_new(:,3);
RAS_Table = table(RAS_labels(2:end), R, A, S);

R = RAS_coords(:,1);
A = RAS_coords(:,2);
S = RAS_coords(:,3);
RAS_Table_Bkp = table(RAS_labels(2:end), R, A, S);

%% Save new coordinates
writetable( RAS_Table_Bkp, fullfile( reconPath, Pt, [Pt '_RAS_Backup.xlsx'] ) )
writetable( RAS_Table, fullfile( reconPath, Pt, [Pt '_RAS.xlsx'] ) )

        
