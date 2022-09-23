%% Template Matlab script to create an BIDS compatible sub-01_ses-01_task-FullExample-01_events.tsv file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% anushkab, 2018

%%
% clear;
% root_dir = ['..' filesep '..'];
% project_label = 'templates';
% sub_id = '01';
% ses_id = '01';
% task_id = 'FullExample';

acquisition = 'ieeg';
% run_id = '01';


%% make an event table and save

for bnk=1:2
    if bnk==1
        Trials=TrialAlignedMat1;
        
        events_tsv_name = fullfile(root_dir, project_label, ...
            ['sub-' ieeg_sub], ...
            ['ses-' ieeg_ses], acquisition, ...
            ['sub-' ieeg_sub ...
            '_ses-' ieeg_ses ...
            '_task-' ieeg_task ...
            'NSP1'...
            '_run-' ieeg_run '_events.tsv']);
    elseif bnk==2 
        Trials=TrialAlignedMat2;
        
        events_tsv_name = fullfile(root_dir, project_label, ...
            ['sub-' ieeg_sub], ...
            ['ses-' ieeg_ses], acquisition, ...
            ['sub-' ieeg_sub ...
            '_ses-' ieeg_ses ...
            '_task-' ieeg_task ...
            'NSP2'...
            '_run-' ieeg_run '_events.tsv']);
    end
    if isempty(Trials)==0
        
        %% CONTAINS a set of REQUIRED and OPTIONAL columns
        % REQUIRED Onset (in seconds) of the event measured from the beginning of
        % the acquisition of the first volume in the corresponding task imaging data file.
        % If any acquired scans have been discarded before forming the imaging data file,
        % ensure that a time of 0 corresponds to the first image stored. In other words
        % negative numbers in onset are allowed.
        onset = Trials(:,3);
        
        % REQUIRED. Duration of the event (measured from onset) in seconds.
        % Must always be either zero or positive. A "duration" value of zero implies
        % that the delta function or event is so short as to be effectively modeled as an impulse.
        duration=NaN*ones(size(Trials,1),1);
        Frequency=NaN*ones(size(Trials,1),1);
        electrical_stimulation_site=cell(size(Trials,1),1);
        % OPTIONAL Primary categorisation of each trial to identify them as instances
        % of the experimental conditions
        electrical_stimulation_trial_type=cell(size(Trials,1),1);
        for TR=1:size(Trials,1)
            if Trials(TR,6)==0
                durStim=0.233;
                pulsetype='Single Pulse';
                freq=0;
            elseif Trials(TR,6)==1
                durStim=1.053;
                pulsetype='Single Pulse';
                freq=0;
            else
                durStim=400;
                pulsetype='Train';
                freq=Trials(TR,6);
            end
            duration(TR) = durStim;
            Frequency(TR)=freq;
            if Trials(TR,9)==1
                electrical_stimulation_trial_type(TR) = {['Bipolar ',pulsetype,' stim, channels: ',ElectrodeLabelsBank1{Trials(TR,8)},'-',ElectrodeLabelsBank1{Trials(TR,10)},...
                    ', current:',num2str(Trials(TR,7)),' mA',', stimulation duration:',num2str(durStim),'ms, frequency:',num2str(freq),'Hz']};
            elseif Trials(TR,9)==2
                electrical_stimulation_trial_type(TR) =  {['Bipolar ',pulsetype,' stim, channels: ',ElectrodeLabelsBank2{Trials(TR,8)},'-',ElectrodeLabelsBank2{Trials(TR,10)},...
                    ', current:',num2str(Trials(TR,7)),' mA',', stimulation duration:',num2str(durStim),'ms, frequency:',num2str(freq),'Hz']};
            end
            electrical_stimulation_site(TR)= {[ElectrodeLabelsBank1{Trials(TR,8)},'-',ElectrodeLabelsBank1{Trials(TR,10)}]};
        end
        
        % OPTIONAL. Response time measured in seconds. A negative response time can be
        % used to represent preemptive responses and n/a denotes a missed response.
        response_time = repmat({'n/a'},size(Trials,1),1);
        
        % OPTIONAL Represents the location of the stimulus file (image, video, sound etc.)
        % presented at the given onset time
        stim_file = repmat({'n/a'},size(Trials,1),1);
        
        % OPTIONAL Hierarchical Event Descriptor (HED) Tag.
        HED = repmat({'n/a'},size(Trials,1),1);
        
        % OPTIONAL Hierarchical Event Descriptor (HED) Tag.
        electrical_stimulation_type= repmat({'biphasic'},size(Trials,1),1);
        electrical_stimulation_StimOnset_sec = Trials(:,2);
        electrical_stimulation_StimOffset_sec = Trials(:,2);
        electrical_stimulation_Frequency_Hz = Frequency;
        electrical_stimulation_Current_mA = Trials(:,7);
        electrical_stimulation_Channel1Stim = Trials(:,8);
        electrical_stimulation_Channel1StimNSP = Trials(:,9);
        electrical_stimulation_Channel2Stim = Trials(:,10);
        electrical_stimulation_Channel2StimNSP = Trials(:,11);
        %% Save table
        t = table(onset, duration, electrical_stimulation_trial_type, response_time, stim_file, HED,...
            electrical_stimulation_StimOnset_sec,electrical_stimulation_StimOffset_sec,electrical_stimulation_Frequency_Hz,electrical_stimulation_Current_mA,...
            electrical_stimulation_site,electrical_stimulation_type,electrical_stimulation_Channel1Stim,...
            electrical_stimulation_Channel1StimNSP,electrical_stimulation_Channel2Stim,electrical_stimulation_Channel2StimNSP);
        
        writetable(t, events_tsv_name, 'FileType', 'text', 'Delimiter', '\t');
        
        %% associated data dictionary
        
        template = struct( ...
            'LongName', TaskName, ...
            'Description', '', ...
            'Levels', [], ...
            'Units', 'n/a', ...
            'TermURL', 'https://github.com/Center-For-Neurotechnology/CereLAB');
        
        dd_json.trial_type = template;
        dd_json.trial_type.Description = 'Stimulation';
        dd_json.trial_type.Levels = struct( ...
            'SinglePulseshort', '0.233ms duration single pulse stim, multiple current amplitudes, multiple locations', ...
            'SinglePulselong', '1.053ms duration single pulse stim, 7mA current amplitudes, multiple locations', ...
            'Train', '400ms duration, multiple frequencies from 10Hz to 200Hz, multiple locations');
        
        dd_json.identifier.LongName = 'Intracranial direct electrical stimulation';
        dd_json.identifier.Description = ['Stimulation was controlled via a custom Cerestim API via MATLAB or a custom C++ code (https://github.com/Center-For-Neurotechnology/CereLAB).',... 
'Waveforms of two different durations were used: .',...
'233 µs duration: 90 µs charge-balanced biphasic symmetrical pulses with an interphase interval of 53 µsec with between 5 and 100 trials with a median of 20 trials per stimulation site across the data set .',...
'and 2) 1053 µs (~1 msec) duration: 500 µs charge-balanced biphasic symmetrical pulses with an interphase interval of 53 µsec with between 10 and 26 trials per stimulation site with a median of 10 trials per site. .',...
'The interval at 53 ms was required as a hardware-limited minimum interval between square pulses with the CereStim stimulator. .',...
'For varying Current Amplitude tests, multiple current amplitudes were applied with the short duration (233 µsec) bipolar stimulation at the following steps: 0.5 mA to 10 mA at 0.5 mA steps with a minimum of 10 trials per stimulation site. ',...
'For train stimuli, each trial delivered a 400 ms train (90 ms charge-balanced biphasic symmetrical pulses with an interphase interval of 53 msec) with varying ',...
'pulse frequencies (10-200 Hz) and amplitudes (0.5-10 mA).']; % Free form description of stimulation parameters,

        
        dd_json.StimulusPresentation.OperatingSystem = 'Windows 10';
        dd_json.StimulusPresentation.SoftwareName = 'CereStim API and CereLAB';
        dd_json.StimulusPresentation.SoftwareRRID = 'xxx_xxx';
        dd_json.StimulusPresentation.SoftwareVersion = 'n/a';
        dd_json.StimulusPresentation.Code = 'https://github.com/Center-For-Neurotechnology/CereLAB';
        
        %% Write JSON
        
        json_options.indent = ' '; % this just makes the json file look prettier
        % when opened in a text editor
        
        jsonSaveDir = fileparts(events_tsv_name);
        if ~isdir(jsonSaveDir)
            fprintf('Warning: directory to save json file does not exist: %s \n', jsonSaveDir);
        end
        
        try
            jsonwrite(strrep(events_tsv_name, '.tsv', '.json'), dd_json, json_options);
        catch
            warning('%s\n%s\n%s\n%s', ...
                'Writing the JSON file seems to have failed.', ...
                'Make sure that the following library is in the matlab/octave path:', ...
                'https://github.com/gllmflndn/JSONio');
        end
    end
end

