%% Template Matlab script to create an BIDS compatible dataset_description.json file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% Writing json files relies on the JSONio library
% https://github.com/gllmflndn/JSONio
% Make sure it is in the matab/octave path
%
% DHermes, 2017

%%
% clear;
% root_dir = ['..' filesep '..'];
% project_label = 'templates';
root_dir='Y:\StimDataBackup\Data_Stimulation\DeIdentifiedDataSet\';
project_label='';
json_label = 'dataset_description';

dataset_description_json_name = fullfile(root_dir, project_label, ...
                                         'dataset_description.json');

%% General fields, shared with MRI BIDS and MEG BIDS:

%% Required fields:

dd_json.Name = 'Intracranial Single Pulse Stimulation'; % name of the dataset

dd_json.BIDSVersion = '1.7.0'; % The version of the BIDS standard that was used

% The interpretation of the dataset. MUST be one of "raw" or "derivative".
% For backwards compatibility, the default value is "raw".
dd_json.DatasetType = 'raw';

%% Recommended fields:

dd_json.License = ['This dataset (EEG and MRI data) is proprietary of the Massachusetts General Hospital and Brigham and Womens Hospital under the PDDL (Open Data Commons Public Domain Dedication and License) which is the License to assign public domain like permissions without giving up the copyright. ']; % what license is this dataset distributed under? The
% use of license name abbreviations is suggested for specifying a license.
% A list of common licenses with suggested abbreviations can be found in appendix III.

dd_json.Authors = {'Sydney S. Cash','Angelique C. Paulk', 'Rina Zelmann', 'Britni Crocker','Alik S. Widge',...
    'Emad N. Eskandar','Ziv M. Williams','R. Mark Richardson','Daniel S. Weisholtz','G. Rees Cosgrove'}; % List of individuals who contributed to the
% creation/curation of the dataset

dd_json.Acknowledgements = 'We would like to thank Giovanni Piantoni, Jean-Baptiste Eichenlaub, Erica Johnson, Gavin Belok, Mia Borzello, Kara Farnes, Dan Soper, Constantin Krempp, Jaquelin Dezha-Peralta, and Pariya Salami for help in data collection. We would like to especially thank the patients for participating in the study. '; % who should be acknowledge in helping to collect the data

dd_json.HowToAcknowledge = ['This project contains the data for the publication Paulk et al, "Local and distant cortical responses to single pulse intracranial stimulation in the human brain are differentially modulated by specific stimulation parameters". ',...
'It contains the raw and preprocessed (bipolarized and epoched) intracranial EEG (iEEG) data files for 52 participants (subjects) with intractable epilepsy implanted with intracranial electrodes for clinically indicated purposes. The data set involves the use of direct electrical stimulation to examine effects of stimulation in the brain. ',...
'Contact: Angelique C. Paulk (apaulk@mgh.harvard.edu) ',...
'If you use this data as a part of any publications, please use the following citation: ',...
'Paulk AC, Zelmann R, Crocker B, Widge AS, Dougherty DD, Eskandar EN, Weisholtz DS, Richardson RM, Cosgrove GR, Williams ZM, Cash SS (2022) Local and distant cortical responses to single pulse intracranial stimulation in the human brain are differentially modulated by specific stimulation parameters. Brain Stimul Available at: https://www.sciencedirect.com/science/article/pii/S1935861X22000456.']; % Instructions how researchers using this
% dataset should acknowledge the original authors. This field can also be used
% to define a publication that should be cited in publications that use the
% dataset.

dd_json.Funding = {'NIH MH086400', 'NIH DA026297', 'NIH EY017658',...
    'NIH MH109722','NIH NS100548','NIH MH111872','NIH NS100548','NIH K24-NS088568',...
    'Tiny Blue Dot Foundation','DE-FG02-97ER25308','DARPA W911NF-14-2-0045'}; % sources of funding (grant numbers)

% List of ethics committee approvals of the research protocols and/or protocol identifiers.
dd_json.EthicsApprovals = {'All patients voluntarily participated after fully informed consent as monitored by the Partners Institutional Review Board covering Brigham and Women’s Hospital (BWH) and Massachusetts General Hospital (MGH). Participants were informed that participation in the stimulation tests would not alter their clinical treatment in any way, and that they could withdraw at any time without jeopardizing their clinical care. '};

% a list of references to
% publication that contain information on the dataset, or links.
dd_json.ReferencesAndLinks = {'Paulk AC, Zelmann R, Crocker B, Widge AS, Dougherty DD, Eskandar EN, Weisholtz DS, Richardson RM, Cosgrove GR, Williams ZM, Cash SS (2022) Local and distant cortical responses to single pulse intracranial stimulation in the human brain are differentially modulated by specific stimulation parameters. Brain Stimul Available at: https://www.sciencedirect.com/science/article/pii/S1935861X22000456. ',...
    'Crocker B, Ostrowski L, Williams ZM, Dougherty DD, Eskandar EN, Widge AS, Chu CJ, Cash SS, Paulk AC (2021) Local and Distant responses to single pulse electrical stimulation reflect different forms of connectivity. Neuroimage 237:118094 Available at: https://www.sciencedirect.com/science/article/pii/S1053811921003712.',...
    'Basu I, Robertson MM, Crocker B, Peled N, Farnes K, Vallejo-Lopez DI, Deng H, Thombs M, Martinez-Rubio C, Cheng JJ, McDonald E, Dougherty DD, Eskandar EN, Widge AS, Paulk AC, Cash SS (2019) Consistent Linear and Non-Linear Responses to Electrical Brain Stimulation Across Individuals and Primate Species. Brain Stimul 12:877–892 Available at: https://doi.org/10.1016/j.brs.2019.03.007. '};

dd_json.DatasetDOI = 'https://doi.org/10.18120/j3k5-1m18'; % the Document Object Identifier of the dataset
% (not the corresponding paper).

%% Write JSON

% this just makes the json file look prettier
% when opened in a text editor
json_options.indent = ' ';

jsonSaveDir = fileparts(dataset_description_json_name);
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist: %s \n', jsonSaveDir);
end

try
    jsonwrite(dataset_description_json_name, dd_json, json_options);
catch
    warning('%s\n%s\n%s\n%s', ...
            'Writing the JSON file seems to have failed.', ...
            'Make sure that the following library is in the matlab/octave path:', ...
            'https://github.com/gllmflndn/JSONio');
end