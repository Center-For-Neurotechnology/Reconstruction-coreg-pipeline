%% Template Matlab script to create an BIDS compatible electrodes.tsv file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% DHermes, 2017
% modified Jaap van der Aar 30.11.18

% clear;

% ScalpElectrodeList={'FP1';'F7';'T3';'T5';'O1';'F3';'C3';'P3';'FP2';'F8';'T4';'T6';'O2';'F4';'P4';'FZ';'CZ';'PZ';'T1';'T2';'CII'};
% EOGElectrodesList={'LOC';'ROC';};
% EMGElectrodesList={'EMG'};
% ECGElectrodesList={'EKG'};
% MiscElectrodeslist={'RESP','chan'};
% TrigElectrodesList={{'SYNC';'TRIGGER';'IMAGE';'34';'DETECT';'SHAM';'SEND STIM';'38';'39';'40';'41';'43';'44';'45';'46'}};
%%

% root_dir = ['Y:\StimDataBackup\Data_Stimulation\MG141\NetworkWake1\FT_Analyzed\UpdatedFT_Analyzed\bidsiEEG2\'];
% project_label = 'SinglePulseStimulation';
% ieeg_sub = '01';
% ieeg_ses = '01';
% ieeg_task = StimInfo;
% ieeg_run = '01';
% 
% mkdir(fullfile(root_dir, project_label,['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg'));
% 
% load(FileLoad,'ElectrodeLabelsBank1','ElectrodeLabelsBank2')
% load(ChannelBank1)
% load(ChannelBank2)
%% make a participants table and save



%% required columns

for bnk=1:2
    if bnk==1
        ChannelList=ElectrodeLabelsBank1;
        %         GoodBadChan1=zeros(size(ChannelList,1),1);
        %         GoodBadChan1(BadChan1)=1;
        
        channels_tsv_name = fullfile(root_dir, project_label, ...
            ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
            ['sub-' ieeg_sub ...
            '_ses-' ieeg_ses ...
            '_task-' ieeg_task ...
            'NSP1'...
            '_run-' ieeg_run '_channels.tsv']);
    elseif bnk==2
        ChannelList=ElectrodeLabelsBank2;
        %         GoodBadChan1=zeros(size(ChannelList,1),1);
        %         GoodBadChan1(BadChan2)=1;
        
        channels_tsv_name = fullfile(root_dir, project_label, ...
            ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
            ['sub-' ieeg_sub ...
            '_ses-' ieeg_ses ...
            '_task-' ieeg_task ...
            'NSP2'...
            '_run-' ieeg_run '_channels.tsv']);
    end
    if isempty(ChannelList)==0
        type=ChannelInformation(ChannelPairs(:,3)==bnk,12);
        name=ChannelInformation(ChannelPairs(:,3)==bnk,11);
        units=ChannelInformation(ChannelPairs(:,3)==bnk,13);
        low_cutoff=ChannelInformation(ChannelPairs(:,3)==bnk,15);
        high_cutoff=ChannelInformation(ChannelPairs(:,3)==bnk,14);
        reference=ChannelInformation(ChannelPairs(:,3)==bnk,16);
        group=ChannelInformation(ChannelPairs(:,3)==bnk,17);
        sampling_frequency=ChannelInformation(ChannelPairs(:,3)==bnk,18);
        status_description=ChannelInformation(ChannelPairs(:,3)==bnk,22);
        for simw=1:length(status_description)
            if isempty(status_description{simw})==1
                status_description{simw}='n/a';
            end
        end
        description=ChannelInformation(ChannelPairs(:,3)==bnk,19);
        status=ChannelInformation(ChannelPairs(:,3)==bnk,21);
        for simw=1:length(status)
            if isempty(status{simw})==1
                status{simw}='n/a';
            end
        end
        notch=ChannelInformation(ChannelPairs(:,3)==bnk,20);
        
        
        %% write
        t = table(name, type, units, low_cutoff, high_cutoff, reference, ...
            group, sampling_frequency, description, notch, status, status_description);
        if bnk==1
            t1=t;
        elseif bnk==2
            t2=t;
        end
        %
        writetable(t, channels_tsv_name, 'FileType', 'text', 'Delimiter', '\t');
    end
end

