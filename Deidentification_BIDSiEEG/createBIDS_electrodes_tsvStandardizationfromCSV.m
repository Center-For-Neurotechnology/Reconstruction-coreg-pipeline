%% Template Matlab script to create an BIDS compatible electrodes.tsv file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% DHermes, 2017
% modified Giulio Castegnaro 201811
% modified further

% Simple bid of code to calculate the surface area of electrode contacts
% r=.4;
% h=2;
% A=2*pi*r*h
%
clear size Size;

ieeg_sub='5o1r';
ieeg_ses = 'postimp';
wd = ['Y:\ReconPipelinePaper\Data\'];

PostopMRIDirectory = fullfile(wd , [PatientName], 'ses-postimp', 'ieeg'); 
RASFile=dir([PostopMRIDirectory,'\*.csv']);
ChInfo=readtable([PostopMRIDirectory,'\',RASFile(1).name]);
ChannelInformation=table2cell(ChInfo);

electrodes_tsv_name = fullfile(wd, ...
    ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
    ['sub-' ieeg_sub ...
    '_ses-' ieeg_ses ...
    '_electrodes.tsv']);

%% make a participants table and save
ElectrodeInfo={};
svc=1;
RAS=[];
for sismn=1:size(ChannelInformation,1)
    if isempty(ChannelInformation{sismn,1})==0 
        RAS(svc,:)=[ChannelInformation{sismn,2} ChannelInformation{sismn,3} ChannelInformation{sismn,4} ChannelInformation{sismn,7}];
        ElectrodeInfo(svc,:)=ChannelInformation(sismn,:);
        svc=1+svc;
    end
end
%% required columns
name = ElectrodeInfo(:,1); % Name of the electrode contact point
x = RAS(:,1); % X position. The positions of the center of each electrode in xyz space.
% Units are in millimeters or pixels and are specified in _*space-<label>_electrode.json.
y = RAS(:,2); % Y position.
z = RAS(:,3); % Z position. If electrodes are in 2D space this should be a column of n/a values.
Size = cell2mat(ElectrodeInfo(:,7)); % Surface area in mm^2 A=2πrh+2πr2

%% recommended columns
material = repmat({'platinum'},size(ElectrodeInfo,1),1); % Material of the electrodes
manufacturer = ElectrodeInfo(:,6); % Optional field to specify the electrode manufacturer
% for each electrode. Can be used if electrodes were manufactured by more than one company.
group = ElectrodeInfo(:,8); % Optional field to specify the group that the electrode is a part of.
% Note that any group specified here should match a group specified in `_channels.tsv`
hemisphere = ElectrodeInfo(:,9); % Optional field to specify the hemisphere in which
% the electrode is placed, one of ["L" or "R"] (use capital).

%% optional columns
type = repmat({'n/a'},size(ElectrodeInfo,1),1); % Optional type of the electrode,for example:cup, ring, clip-on, wire, needle, ...
impedance = repmat({'n/a'},size(ElectrodeInfo,1),1); % Impedance of the electrode in kOhm.

%% write a tsv file and a csv file just for later use
size=Size;
t = table(name, x, y, z, size, material, manufacturer, group, hemisphere, type, impedance);

writetable(t, electrodes_tsv_name, 'FileType', 'text', 'Delimiter', '\t');

writetable(t,[electrodes_tsv_name(1:end-4),'.csv'])
clear size Size;
