%% Load the inputs for the snap2surf file
Pt = 'id##';
ElectrodeRASPath = fullfile( '/home/bourbon-the-huckster/Documents/Recons', Pt, [Pt '_RAS.xlsx']); 
[num,txt,all]=xlsread(ElectrodeRASPath);

%Number of RAS Coords
RAS_coords=[];
RAS_coords=num(1:length(num),:);

RAS_labels={};
RAS_labels=txt(1:length(txt),1);

%% Separate each of the grids and strips into different variables
%G = RAS_coords(1:56,:);%Lateral Temporal Grid
gr1 = RAS_coords(1:8,:);
gr2 = RAS_coords(9:16,:);
gr3 = RAS_coords(17:24,:);
gr4 = RAS_coords(25:32,:);

pgr1 = RAS_coords(33:37,:);
pgr2 = RAS_coords(38:42,:);
pgr3 = RAS_coords(43:47,:);
pgr4 = RAS_coords(48:52,:);

REF = RAS_coords(79:82,:);

%Depths -- Will have to be adjusted manually below if you wish to move them
depths = RAS_coords(83:98,:);


%% Code for adjusting the locations of the grids and strips
surface_path = fullfile( '/home/bourbon-the-huckster/Documents/Recons' , Pt, [Pt '_SurferOutput/surf'] );
% side = 'r';
% surftype = 'pial-outer-smoothed';
[surf.vert surf.tri]=read_surf([surface_path '/rh.pial-outer-smoothed']);

[gr1_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(gr1,surf,nchoosek(1:32,2));
[gr2_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(gr2,surf);
[gr3_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(gr3,surf);
[gr4_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(gr4,surf);

[pgr1_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(pgr1,surf);%,nchoosek(1:32,2)
[pgr2_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(pgr2,surf);
[pgr3_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(pgr3,surf);
[pgr4_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(pgr4,surf);

[REF_coor_snapped,GR_electrode_pairouts] = snap2dural_energy_customized(REF,surf);

%% Compile new coordinates
RAS_coords_new = [gr1_coor_snapped;gr2_coor_snapped;gr3_coor_snapped;gr4_coor_snapped;pgr1_coor_snapped;...
    pgr2_coor_snapped; pgr3_coor_snapped; pgr4_coor_snapped; ATS_coor_snapped; PTS_coor_snapped;...
    AIS_coor_snapped; CIS_coor_snapped; PIS_coor_snapped; REF_coor_snapped; depths ];
RAS_as_cell = num2cell(RAS_coords_new);
R = RAS_coords_new(:,1);
A = RAS_coords_new(:,2);
S = RAS_coords_new(:,3);
RAS_Table = table(RAS_labels(2:end), R, A, S);


%% Save new coordinates
writetable( RAS_Table, fullfile( '/home/bourbon-the-huckster/Documents/Recons', Pt, [Pt '_RAS_Snapped.xlsx'] ) )

        
