%% General fields, shared with MRI BIDS and MEG BIDS:
                    ieeg_json_name = fullfile(root_dir,project_label, ['sub-' ieeg_sub], ...
                        ['ses-' ieeg_ses], ...
                        'ieeg', ...
                        ['sub-' ieeg_sub ...
                        '_ses-' ieeg_ses ...
                        '_task-' ieeg_task ...
                        '_run-' ieeg_run '_ieeg.json']);
%% Required fields:
ieeg_json.TaskName = TaskName; % Name of the task (for resting state use the rest
% prefix). No two tasks should have the same name. Task label is derived
% from this field by removing all non alphanumeric ([a-zA-Z0-9]) characters.
% Note this does not have to be a "behavioral task" that subjects perform, but can reflect some
% information about the conditions present when the data was acquired (for example: "rest" or "sleep").

ieeg_json.SamplingFrequency = FS1; % Sampling frequency (in Hz) of all the iEEG channels
% in the recording (for example: 2400). All other channels should have frequency specified
% as well in the channels.tsv file.

ieeg_json.PowerLineFrequency = 60; % Frequency (in Hz) of the power grid where the
% iEEG recording was done (for instance: 50 or 60).

ieeg_json.SoftwareFilters = 'n/a'; % List of temporal software filters applied or
% ideally key:value pairs of pre-applied filters and their parameter values.
% (n/a if none). for example: "{'HighPass': {'HalfAmplitudeCutOffHz': 1, 'RollOff: '6dB/Octave'}}".

%% Recommended fields:

HardwareFilters.HighpassFilter.CutoffFrequency = 0.3; % Contains the high pass hardware filter
% cut off frequency.

HardwareFilters.LowpassFilter.CutoffFrequency = 7500; % Contains the low pass hardware filter
% cut off frequency.

ieeg_json.HardwareFilters = HardwareFilters; % Cutoff frequencies of high and low pass filter
% are stored in this variable automatically. No further input necessary.

ieeg_json.Manufacturer = 'Blackrock Microsystems'; % Manufacturer of the amplifier system (for example: "TDT, blackrock")

ieeg_json.ManufacturersModelName = 'Neural Signal Processor (Cerebus), FrontEnd amplifiers, CereStim stimulator'; % Manufacturer's designation of the
% iEEG amplifier model (for example: "TDT").

ieeg_json.TaskDescription = 'Stimulation was delivered with a CereStim stimulator (Blackrock Microsystems, Salt Lake City, UT) using single pulse electrical stimulation (SPES). Current injection and return paths used neighboring contacts in a bipolar configuration '; % Longer description of the task.

ieeg_json.Instructions = 'A trained electroencephalographer examined ongoing recordings for epileptiform activity and asked participants if they experienced any sensations. The participants were awake and were aware that they were being stimulated but were blind to the stimulation timing and parameters.'; % Text of the instructions given to participants
% before the recording. This is especially important in context of resting
% state and distinguishing between eyes open and eyes closed paradigms.

ieeg_json.CogAtlasID = 'https://www.sciencedirect.com/science/article/pii/S1935861X22000456'; % URL of the corresponding Cognitive Atlas Task term

ieeg_json.CogPOID = 'https://www.sciencedirect.com/science/article/pii/S1935861X22000456'; % URL of the corresponding CogPO term

ieeg_json.InstitutionName = 'Massachusetts General Hospital'; % The name of the institution in charge of
% the equipment that produced the composite instances.

ieeg_json.InstitutionAddress = '55 Fruit St., Boston, MA 02114, USA'; % The address of the institution in charge
% of the equipment that produced the composite instances.

ieeg_json.DeviceSerialNumber = 'n/a'; % The serial number of the equipment that
% produced the composite instances. A pseudonym can also be used to prevent
% the equipment from being identifiable, as long as each pseudonym is unique
% within the dataset.

ieeg_json.ECOGChannelCount = [(ECogChanNum)]; % Number of iEEG surface channels included in the recording (for example: 120)
ieeg_json.SEEGChannelCount = [(sEEGChanNum)]; % Number of iEEG depth channels included in the recording (for example: 8)
ieeg_json.EEGChannelCount = [(EEGChanNum)]; % Number of scalp EEG channels recorded simultaneously (for example: 21)
ieeg_json.EOGChannelCount = [(EOGChanNum)]; % Number of EOG channels
ieeg_json.ECGChannelCount = [(ECGChanNum)]; % Number of ECG channels
ieeg_json.EMGChannelCount = [(EMGChanNum)]; % Number of EMG channels
ieeg_json.MiscChannelCount = [(MiscChannelCount)]; % Number of miscellaneous analog channels for auxiliary signals
ieeg_json.TriggerChannelCount = [round(TriggerChannelCount)]; % Number of channels for digital (TTL bit level) triggers.
ieeg_json.RecordingDuration = [(RecordingDuration)/FS1]; % Length of the recording in seconds (for example: 3600)
ieeg_json.RecordingType = RecordingType; % Defines whether the recording is "continuous" or "epoched"; this latter
% limited to time windows about events of interest (for example: stimulus presentations, subject responses etc.)
if strmatch(RecordingType,'epoched')==1
ieeg_json.EpochLength = [num2str(EpochLength)]; % Duration of individual epochs in seconds (for example: 1).
end
% If recording was continuous, leave out the field.
ieeg_json.SubjectArtefactDescription = 'n/a'; % Freeform description of the observed
% subject artefact and its possible cause (for example: door open, nurse walked into room at 2 min,
% "Vagus Nerve Stimulator", non-removable implant, seizure at 10 min).
% If this field is left empty, it will be interpreted as absence of artifacts.

% ieeg_json.SoftwareVersions = ''; % Manufacturer's designation of the acquisition software.

%% Specific iEEG fields:

% If mixed types of references, manufacturers or electrodes are used, please
% specify in the corresponding table in the _electrodes.tsv file

%% Required fields:

ieeg_json.iEEGReference = iEEGReference; % General description of the reference scheme used and
% (when applicable) of location of the reference electrode in the raw recordings
% (for example: "left mastoid", "bipolar", "T01" for electrode with name T01, "intracranial electrode
% on top of a grid, not included with data", "upside down electrode"). If different channels have
% a different reference, this field should have a general description and the channel specific
% reference should be defined in the _channels.tsv file.

%% Recommended fields:

ieeg_json.ElectrodeManufacturer = ElectrodeManufacturer; % can be used if all electrodes are of the same manufacturer
% (for example: AD-TECH, DIXI). If electrodes of different manufacturers are used, please use the corresponding
% table in the _electrodes.tsv file
% ieeg_json.ElectrodeManufacturersModelName = ''; % Specify model name. If different electrode types are used,
% please use the corresponding table in the _electrodes.tsv file
ieeg_json.iEEGGround = iEEGGround; % Description of the location of the ground electrode
% ("placed on right mastoid (M2)").
ieeg_json.iEEGPlacementScheme = iEEGPlacementScheme; % Freeform description of the placement of the iEEG electrodes.
% Left/right/bilateral/depth/surface(for example: "left frontal grid and bilateral hippocampal depth" or
% "surface strip and STN depth" or "clinical indication bitemporal, bilateral temporal strips and left grid").
ieeg_json.iEEGElectrodeGroups = GroupList; % Field to describe the way electrodes are grouped
% into strips, grids or depth probes for example: {'grid1': "10x8 grid on left temporal pole",
% 'strip2': "1x8 electrode strip on xxx"}.

%% Optional fields:

ieeg_json.ElectricalStimulation = true; % Boolean field to specify if electrical stimulation
% was done during the recording (options are "true" or "false"). Parameters for event-like
% stimulation should be specified in the _events.tsv file (see example underneath).

ieeg_json.ElectricalStimulationParameters = ['Stimulation was controlled via a custom Cerestim API via MATLAB or a custom C++ code (https://github.com/Center-For-Neurotechnology/CereLAB).',... 
'Waveforms of two different durations were used: .',...
'233 µs duration: 90 µs charge-balanced biphasic symmetrical pulses with an interphase interval of 53 µsec with between 5 and 100 trials with a median of 20 trials per stimulation site across the data set .',...
'and 2) 1053 µs (~1 msec) duration: 500 µs charge-balanced biphasic symmetrical pulses with an interphase interval of 53 µsec with between 10 and 26 trials per stimulation site with a median of 10 trials per site. .',...
'The interval at 53 ms was required as a hardware-limited minimum interval between square pulses with the CereStim stimulator. .',...
'For varying Current Amplitude tests, multiple current amplitudes were applied with the short duration (233 µsec) bipolar stimulation at the following steps: 0.5 mA to 10 mA at 0.5 mA steps with a minimum of 10 trials per stimulation site. ']; % Free form description of stimulation parameters,
% such as frequency, shape etc. Specific onsets can be specified in the _events.tsv file.
% Specific shapes can be described here in freeform text.

%% write
jsonSaveDir = fileparts(ieeg_json_name);
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist, create: %s \n', jsonSaveDir);
end

json_options.indent = ' '; % this just makes the json file look prettier
% when opened in a text editor

try
    jsonwrite(ieeg_json_name, ieeg_json, json_options);
catch
    warning('%s\n%s\n%s\n%s', ...
            'Writing the JSON file seems to have failed.', ...
            'Make sure that the following library is in the matlab/octave path:', ...
            'https://github.com/gllmflndn/JSONio');
end