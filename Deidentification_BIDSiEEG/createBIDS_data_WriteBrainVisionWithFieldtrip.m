%
% %%%%% save data as BrainVision BIDS %%%%%
%
% Sample script that calls Fieldtrip functions to write a Brainvision dataset
% Added fields are examples, read these from the raw data
%
% Fieldtrip has to be in the path!
%
% D. Hermes 2018

% provide necessary labels and path to save the data:
% dataRootPath = 'Y:\StimDataBackup\Data_Stimulation\MG141\NetworkWake1\FT_Analyzed\UpdatedFT_Analyzed\bidsiEEG2\';
% sub_label = '01';
% ses_label = '01';
% task_label = 'networkWake1';
% run_label = '01';
%
% % name to save data:


for bnk=1:2
    
    if bnk==1
        DataFeed=DATAlfp1;
        
        ieeg_name_save = fullfile(root_dir,project_label, ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
            ['sub-' ieeg_sub...
            '_ses-' ieeg_ses...
            '_task-' ieeg_task...
            'NSP1'...
            '_run-' ieeg_run...
            '_ieeg']);
    elseif bnk==2
        DataFeed=DATAlfp2;
        ieeg_name_save = fullfile(root_dir,project_label, ['sub-' ieeg_sub], ['ses-' ieeg_ses], 'ieeg', ...
            ['sub-' ieeg_sub...
            '_ses-' ieeg_ses...
            '_task-' ieeg_task...
            'NSP2'...
            '_run-' ieeg_run ...
            '_ieeg']);
    end
    
    
    if isempty(DataFeed)==0
    %%%% assign header fields:
    
    % sampling frequency
    dataStruct.hdr.Fs = FS1;
    dataStruct.hdr.nChans = size(DataFeed,1);
    
    % number of channels
    dataStruct.hdr.label = cell(dataStruct.hdr.nChans, 1);
    for kk = 1:dataStruct.hdr.nChans
        dataStruct.hdr.label{kk} = ['iEEG' int2str(kk)];
    end
    
    % number of samples
    dataStruct.hdr.nSamples = size(DataFeed,2);
    
    % ?
    dataStruct.hdr.nSamplesPre = 0;
    
    % 1 trial for continuous data
    dataStruct.hdr.nTrials = 1;
    
    % channels type, see BIDS list of types
    dataStruct.hdr.chantype = cell(dataStruct.hdr.nChans, 1);
    for kk = 1:124
        dataStruct.hdr.chantype{kk} = ['sEEG'];
    end
    % for kk = 125:128
    %     dataStruct.hdr.chantype{kk} = ['scalp'];
    % end
    % for kk = 129:143
    %     dataStruct.hdr.chantype{kk} = ['triggers'];
    % end
    
    % I still don't how to indicate the mu in BIDS, now using letter u
    dataStruct.hdr.chanunit = repmat({'uV'}, size(dataStruct.hdr.chantype, 1), 1);
    
    % labels again, same as before
    dataStruct.label = dataStruct.hdr.label;
    
    % time vector
    dataStruct.time{1} = [1:dataStruct.hdr.nSamples] / dataStruct.hdr.Fs;
    
    % put the data matrix here: electrodes x samples
    dataStruct.trial{1} = DataFeed;
    
    % sampling freq again, same as before
    dataStruct.fsample = dataStruct.hdr.Fs;
    
    % ?
    dataStruct.sampleinfo = [1 dataStruct.hdr.nSamples];
    
    % now fetch a header
    hdr_data = ft_fetch_header(dataStruct);
    
    % save the data
    ft_write_data(ieeg_name_save, dataStruct.trial{1}, 'header', hdr_data, 'dataformat', 'brainvision_eeg');
    end
end
