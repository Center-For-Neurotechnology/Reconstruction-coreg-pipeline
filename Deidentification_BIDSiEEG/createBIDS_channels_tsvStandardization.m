%% Template Matlab script to create an BIDS compatible electrodes.tsv file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% DHermes, 2017
% modified Jaap van der Aar 30.11.18

% clear;

ScalpElectrodeList={'FP1';'F7';'T3';'T7';'T5';'O1';'F3';'C3';'P3';'FP2';'F8';'T4';'T6';'O2';'F4';'P4';'FZ';'CZ';'PZ';'T1';'T2';'CII';'C4';'FP__01';'FP__02';'O3';'C2';'O2 Sa';'FP1_r';'T4_re';'O2';'O2Sat'};
EOGElectrodesList={'LOC';'ROC';};
EMGElectrodesList={'EMG';'EMG1';'EMG2';};
ECGElectrodesList={'EKG';'EKG  ';'Pulse';'EKG1';};
MiscElectrodeslist={'RESP';'chan';'empty';'broke';'104';'105';'106';'107';'108';'109';'110';'111';'112';'113';'114';'115';'116';'117';'118';'119';'120';'121';'122';'123';'124';'125';'126';'127';'128';'BLANK';'BLA1';};
TrigElectrodesList={'SYNC';'TRIGGER';'IMAGE';'34';'DETECT';'SHAM';'SEND STIM';'38';'39';'40';'41';'43';'44';'45';'46';'ainp01';'ainp02';'ainp03';'ainp04';'ainp05';'ainp06';'ainp07';'ainp08';'ainp09';'ainp10';'ainp11';'ainp12';'ainp13';'ainp14';'ainp15';'ainp16';...
    'ainp_1';'ainp_2';'ainp_3';'ainp_4';'ainp_5';'ainp_6';'ainp_7';'ainp_8';'ainp_9';'ainp1';'ainp2';'ainp3';'ainp4';'ainp5';'ainp6';'ainp7';'ainp8';'ainp9';'AUDIO';'TRIGG';'DETEC';...
    'SEND ';'ainp0';'EVT';'trig';'RANDOM';'RANDO';'EVENT';'sync';'REF';'A1';'Trigg';'A1';'A2';'Stim';'SEND';'TRIGGE';'SEND S';'imageOnset';'ainp17';'ainp18';'ainp19';'ainp20';'ainp21';'ainp22';'ainp23';'ainp24';'ainp25';'ainp26';'ainp27';'ainp28';...
    'REF1';'REF3';'REF_01';'REF_03';'BLA0N';'Stim 1';'Stim 2';'Stim1';'Stim2';'cerestimOUT';'alphatrig';'REF_02';'REF_03'};
%%

% root_dir = ['Y:\StimDataBackup\Data_Stimulation\MG141\NetworkWake1\FT_Analyzed\UpdatedFT_Analyzed\bidsiEEG2\'];
% project_label = 'SinglePulseStimulation';
% ieeg_sub = '01';
% ieeg_ses = '01';
% ieeg_task = StimInfo;
% ieeg_run = '01';

% mkdir(fullfile(root_dir, project_label,['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg'));

% load(FileLoad,'ElectrodeLabelsBank1','ElectrodeLabelsBank2')
% load(ChannelBank1)
% load(ChannelBank2)
%% make a participants table and save

BIPOL=0;
PerBankInfo={};
%% required columns

for bnk=1:2
    if bnk==1
        ChannelList=ElectrodeLabelsBank1;
        GoodBadChan=zeros(size(ChannelList,1),1);
        GoodBadChan(BadChan1)=1;
        
%         channels_tsv_name = fullfile(root_dir, project_label, ...
%             ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
%             ['sub-' ieeg_sub ...
%             '_ses-' ieeg_ses ...
%             '_task-' ieeg_task ...
%             '_run-' ieeg_run '_channelsNSP1.tsv']);
    elseif bnk==2
        ChannelList=ElectrodeLabelsBank2;
        GoodBadChan=zeros(size(ChannelList,1),1);
        if isempty(GoodBadChan)==0
        GoodBadChan(BadChan2)=1;
        end
        
%         channels_tsv_name = fullfile(root_dir, project_label, ...
%             ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
%             ['sub-' ieeg_sub ...
%             '_ses-' ieeg_ses ...
%             '_task-' ieeg_task ...
%             '_run-' ieeg_run '_channelsNSP2.tsv']);
    end
    
    type=cell(size(ChannelList,1),1);name=cell(size(ChannelList,1),1);units=cell(size(ChannelList,1),1);
    low_cutoff=cell(size(ChannelList,1),1);high_cutoff=cell(size(ChannelList,1),1);
    reference=cell(size(ChannelList,1),1);
    group=cell(size(ChannelList,1),1);sampling_frequency=cell(size(ChannelList,1),1);
    status_description=cell(size(ChannelList,1),1);description=cell(size(ChannelList,1),1);
    status=cell(size(ChannelList,1),1);notch=cell(size(ChannelList,1),1);
    for chanStep=1:size(ChannelList,1)
        name{chanStep} = ChannelList{chanStep}; % Label of the channel, only contains letters and numbers. The label must
        % correspond to _electrodes.tsv name and all ieeg type channels are required to have \
        % a position. The reference channel name MUST be provided in the reference column
CS=ChannelList{chanStep};
            CS=deblank(ChannelList{chanStep});
            
        for sE=1:length(ScalpElectrodeList)
            if strcmpi(CS,ScalpElectrodeList{sE})
                type{chanStep} = 'EEG'; % Type of channel, see below for adequate keywords in this field
            end
        end
        
        for sE=1:length(ECGElectrodesList)
            if contains(CS,ECGElectrodesList{sE},'IgnoreCase',true)
                type{chanStep} = 'ECG'; % Type of channel, see below for adequate keywords in this field
            end
        end
        
        for sE=1:length(EOGElectrodesList)
            if contains(CS,EOGElectrodesList{sE},'IgnoreCase',true)
                type{chanStep} = 'EOG'; % Type of channel, see below for adequate keywords in this field
            end
        end
        
        for sE=1:length(MiscElectrodeslist)
            if contains(CS,MiscElectrodeslist{sE},'IgnoreCase',true)
                type{chanStep} = 'MISC'; % Type of channel, see below for adequate keywords in this field
            end
        end
        
        for sE=1:length(TrigElectrodesList)
            if strcmpi(CS,TrigElectrodesList{sE})
                type{chanStep} = 'TRIG'; % Type of channel, see below for adequate keywords in this field
            end
        end
        
        for sE=1:length(EMGElectrodesList)
            if contains(CS,EMGElectrodesList{sE},'IgnoreCase',true)
                type{chanStep} = 'EMG'; % Type of channel, see below for adequate keywords in this field
            end
        end
        if isempty(type{chanStep})==1
            type{chanStep} = 'SEEG';
        end
        % Must be one of: "MEGMAG", "MEGGRADAXIAL", "MEGGRADPLANAR", "MEGREFMAG", "MEGREFGRADAXIAL",...
        %     "MEGREFGRADPLANAR", "MEGOTHER", "EEG", "ECOG", "SEEG", "DBS", "VEOG", "HEOG", "EOG", "ECG",...
        %     "EMG", "TRIG", "AUDIO", "PD", "EYEGAZE", "PUPIL", "MISC", "SYSCLOCK", "ADC", "DAC", "HLU", "FITERR", "OTHER".
        
        units{chanStep} = 'uV'; % Physical unit of the value represented in this channel, for example: V for Volt,
        % specified according to the SI unit symbol and possibly prefix symbol (for example: mV, ?V),
        % see the BIDS spec (section 15 Appendix V: Units) for guidelines for Units and Prefixes.
        
        low_cutoff{chanStep} = '1000'; % Frequencies used for the low pass filter applied to the
        % channel in Hz. If no low pass filter was applied, use n/a. Note that
        % anti-alias is a low pass filter, specify its frequencies here if applicable.
        
        high_cutoff{chanStep} = 'n/a'; % Frequencies used for the high pass filter applied to
        % the channel in Hz. If no high pass filter applied, use n/a.
        
        %% recommended columns:
        
        reference{chanStep} = 'mastoid'; % Specification of the reference (for example: "mastoid", "ElectrodeName01",
        % "intracranial", "CAR", "other", "n/a"). If the channel is not an electrode channel
        % (for example: a microphone channel) use `n/a`.
        CL=ChannelList{chanStep};
        CL(CL=='_')=[];
        if strmatch(type{chanStep},'SEEG')==1
            NumOn=[];
            for si=1:length(CL)
                if isempty(str2num(CL(si)))==0
                    if isreal(str2num(CL(si)))==1
                        NumOn=[NumOn;si];
                    end
                end
            end
            group{chanStep}=CL(1:NumOn(1)-1);
        elseif strmatch(type{chanStep},'ECOG')==1
            group{chanStep}=CL(1:NumOn(1)-1);
        else
            group{chanStep}=CL;
        end
        
        %     group{chanStep} = ''; % Which group of channels (grid/strip/probe) this channel belongs to.
        % One group has one wire and noise can be shared. This can be a name or number.
        % Note that any groups specified in `_electrodes.tsv` must match those present here.
        
        %% optional columns
        
        sampling_frequency{chanStep} = '2000'; % Sampling rate of the channel in Hz.
        
        if BIPOL==0
            description{chanStep} = 'referential/unipolar'; % Brief free-text description of the channel, or other information of
        elseif BIPOL==1
            description{chanStep} = 'bipolar'; % Brief free-text description of the channel, or other information of
        end
        % interest (for example: position (for example: 'left lateral temporal surface', 'unipolar/bipolar', etc.)).
        
        notch{chanStep} = 'n/a'; % Frequencies used for the notch filter applied to the channel,
        % in Hz. If no notch filter applied, use n/a.
        
        if isempty(strmatch(type{chanStep},'SEEG'))==0  || isempty(strmatch(type{chanStep},'ECOG'))==0
            if GoodBadChan(chanStep)==0
                status{chanStep} = 'good';
            elseif GoodBadChan(chanStep)==1
                status{chanStep} = 'bad';
            end
        else
            status{chanStep} = 'n/a';
        end
        % Data quality observed on the channel (good/bad). A channel is considered bad
        % if its data quality is compromised by excessive noise. Description of noise type SHOULD be
        % provided in [status_description].
        
        status_description{chanStep} = 'n/a'; % Freeform text description of noise or artifact affecting data
        % quality on the channel. It is meant to explain why the channel was declared bad in [status].
    end
    
    %% write
    t = table(name, type, units, low_cutoff, high_cutoff, reference, ...
        group, sampling_frequency, description, notch, status, status_description);
    if bnk==1
        t1=t;
    elseif bnk==2
       t2=t; 
    end
    PerBankInfo{bnk}=table2cell(t);
    %
%     writetable(t, channels_tsv_name, 'FileType', 'text', 'Delimiter', '\t');
end


