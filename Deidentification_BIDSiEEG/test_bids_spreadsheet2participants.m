% requires:
% - moxunit: https://github.com/MOxUnit/MOxUnit
% - bids-matlab

% addpath(fullfile(pwd, '..'));
pwd='Y:\StimDataBackup\Data_Stimulation\DeIdentifiedDataSet\';
input_file = fullfile(pwd, 'spreadsheet_to_convert.xlsx');
files_out = bids_spreadsheet2participants(input_file, 'ignore', 'comment', 'export', pwd);

json_expected = jsonread(fullfile(pwd, 'data', 'participants.json'));
json_actual = jsonread(fullfile(pwd, 'participants.json'));
% uses moxunit assertEqual
assertEqual(json_actual, json_expected);

tsv_expected  = bids.util.tsvread(fullfile(pwd, 'data', 'participants.tsv'));
tsv_actual  = bids.util.tsvread(fullfile(pwd, 'participants.tsv'));
assertEqual(tsv_actual, tsv_expected);

% A=readtable(input_file,'Sheet','Data')

% sls=7*30*(randn(52,1))

dA=[];
for v = 1:52
    valAdd=[-4 -3 -2 -1 1 2 3 4];
    vk=randperm(6);
    dA(v)=valAdd(vk(1));
end