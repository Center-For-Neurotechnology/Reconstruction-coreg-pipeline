%% iEEG BIDS formatting for converting the stimulation data sets into deidentified
% formatting. This loop takes information from a spreadsheet
clear;
addpath(genpath('X:\Projects\Lab_Materials\Analysis_Tools_and_Software\fieldtrip-20220202'))
addpath(genpath('F:\Dropbox (Personal)\ACPProjects\MechanismsUnderlyingStimulation\GreyWhiteMatter_StimPaper\code'))
addpath(genpath('F:\Dropbox (Personal)\ACPProjects\MechanismsUnderlyingStimulation\GreyWhiteMatter_StimPaper\code\JSONio-main'))

Bidsfolder='Y:\StimDataBackup\Data_Stimulation\DeIdentifiedDataSet\';
PtInfo=readtable(['F:\Dropbox (Personal)\ACPProjects\MechanismsUnderlyingStimulation\Data_Information_Spreadsheets\HPStimFileList_2.xlsx'],'Sheet','MasterList');
Patients=table2cell(PtInfo(1:end,1));
List=table2cell(PtInfo(1:end,3));

PtInfo2=readtable(['F:\Dropbox (Personal)\ACPProjects\MechanismsUnderlyingStimulation\Data_Information_Spreadsheets\HPStimFileList_2.xlsx'],'Sheet','Sheet1');
Patients2=table2cell(PtInfo2(1:end,1));
Identifiers=table2cell(PtInfo2(1:end,3));
Relabel=table2cell(PtInfo2(1:end,4));

Maindir='Y:\StimDataBackup\Data_Stimulation\';
NameTitle={'Network','Amplitude','Stimulation_Matrix','Anesth','Single'};

czv=1;
lsn=2;
SummaryArray=repmat({'NaN'},78,17);
MissingFile=[];MissingFileLabel={};MissingFileSubcort=[];MissingFileLabelSubcort={};NotMissingFileSubcort=[];NotMissingFileLabelSubcort={};
MatchedIDs={};
clow=1;
for DI=75
    if isempty(Identifiers{DI})==0
        czv=1;
        for NT=1:5
            PatientName=Patients{DI};
            PatientNameList=List{DI};
            PatientNameParc=PatientName;
            %
            LoadDir=[Maindir,Patients{DI},'\'];
            
            taskDir=dir([LoadDir,NameTitle{NT},'*']);
            for Tsk=1:length(taskDir)
                TASK=taskDir(Tsk).name;
                
                Taskload=dir([Maindir,PatientName,'\',TASK,'*']);
                
                for session=1:length(Taskload)
                    SessFolder=[Maindir,'\',PatientName,'\',Taskload(session).name,'\'];
                    if isempty(dir([SessFolder,'Cleaned*']))==0
                        MRIDir=dir([Maindir,PatientName,'\participantHP*','_anon.nii*']); % Preop MRI
                        PostopMRIDir=dir([Maindir,PatientName,'\participantHP*','_anonPostOPScan.nii*']); %PostOp MRI
                        FileLoad=dir([SessFolder,'\FT_Analyzed\UpdatedFT_Analyzed\*AlignedToTrialsall','*']); % Raw data saved in Matlab, could be replaced by looking for the .nsp files to open
                        TrialLoad=dir([SessFolder,'\Cleaned','*']); % this had my organized trial information
                        ChannelBank1=dir([SessFolder,'\FT_Analyzed\UpdatedFT_Analyzed\Channel_Selection*','Bank1*']);  % this had my bad channel information, NSP1
                        ChannelBank2=dir([SessFolder,'\FT_Analyzed\UpdatedFT_Analyzed\Channel_Selection*','Bank2*']); % this had my bad channel information, NSP2
                        ParcellationsMapping=dir([SessFolder,'\FT_Analyzed\UpdatedFT_Analyzed\Parcell*','SingleChannel*']); % This is the single channel RAS mapping with the channel information saved from the '' files
                        if isempty(FileLoad)==0 &&  isempty(TrialLoad)==0 &&  isempty(ParcellationsMapping)==0
                            
                            StimInfo=TASK;
                            RawBipolar='RAW';
                            
                            load([SessFolder,'\FT_Analyzed\UpdatedFT_Analyzed\',FileLoad(1).name]) %Loading needed files
                            load([SessFolder,'\FT_Analyzed\UpdatedFT_Analyzed\',ParcellationsMapping(1).name])%Loading needed files
                            load([SessFolder,'\',TrialLoad(1).name],'TrialAlignedMat1','TrialAlignedMat2')%Loading trial information
                            %                     if NT>1
                            %                         pause
                            %                     end
                            
                            %% This section was to fix some few files with weird trials or the wrong column order.
                            %                             I can fix those specific examples later once this needs to be run 'for real'
                            if NT==3
                                TrialAlignedMat1=TrialAlignedMat1(:,[1:5 7 6 8:12]);
                                TrialAlignedMat1(TrialAlignedMat1(:,6)==1)=0;
                                if isempty(TrialAlignedMat2)==0
                                    TrialAlignedMat2=TrialAlignedMat2(:,[1:5 7 6 8:12]);
                                    TrialAlignedMat2(TrialAlignedMat2(:,6)==1)=0;
                                end
                            end
                            if NT==5
                                TrialAlignedMat1=TrialAlignedMat1(:,[1:5 7 6 8:12]);
                                TrialAlignedMat1(TrialAlignedMat1(:,6)==1)=0;
                                if isempty(TrialAlignedMat2)==0
                                    TrialAlignedMat2=TrialAlignedMat2(:,[1:5 7 6 8:12]);
                                    TrialAlignedMat2(TrialAlignedMat2(:,6)==1)=0;
                                end
                            end
                            
                            
                            if NT==5 && DI==36 && Tsk==2
                                TrialAlignedMat1=[TrialAlignedMat1(:,[1:5 7:12]) NaN*ones(size(TrialAlignedMat1,1),1)];
                                TrialAlignedMat1(TrialAlignedMat1(:,6)==1)=0;
                                if isempty(TrialAlignedMat2)==0
                                    TrialAlignedMat2=[TrialAlignedMat2(:,[1:5 7 6 8:12]) NaN*ones(size(TrialAlignedMat1,1),1)];
                                    TrialAlignedMat2(TrialAlignedMat2(:,6)==1)=0;
                                end
                            end
                            %% This section generates the new BIDS file format
                            %                             organization including the imaging folders
                            
                            TaskName=(Taskload(session).name);
                            
                            root_dir = [Bidsfolder,'\'];
                            ieeg_project = 'SinglePulseStimulation';
                            ieeg_sub = Relabel{DI}; %Essential to name everything correctly.
                            ieeg_ses = 'postimp';
                            taskName=[TaskName,'',RawBipolar];
                            taskName(taskName=='-')=[]; %underscores and dashes will break the BIDS system so all task names must have those removed.
                            taskName(taskName=='_')=[]; %underscores and dashes will break the BIDS system so all task names must have those removed.
                            ieeg_task = taskName;
                            if length(num2str(czv))==1
                                ieeg_run = ['0',num2str(czv)];
                            else
                                ieeg_run = [num2str(czv)];
                            end
                            project_label='';
                            
                            %% The variables below are important for different .json coordinate files, etc. 
                            FS1=2000; %Sample rate, change if this depends on the recording
                            chType1=(ChannelInformation(:,12));
                            grouptype1=(ChannelInformation(:,17));
                            
                            ECogChanNum=sum(count(chType1,'ECOG'));
                            sEEGChanNum=sum(count(chType1,'SEEG'));
                            eegonly=count(chType1,'EEG');
                            EEGChanNum=sum(count(chType1,'SEEG')==0 & eegonly==1);
                            EOGChanNum=sum(count(chType1,'EOG'));
                            ECGChanNum=sum(count(chType1,'EKG'));
                            EMGChanNum=sum(count(chType1,'EMG'));
                            MiscChannelCount=sum(count(chType1,'MISC'));
                            TriggerChannelCount=sum(count(chType1,'TRIG'));
                            
                            RecordingDuration=size(DATAlfp1,2);
                            RecordingType='continuous';
                            EpochLength=[];
                            %                     iEEGReference='mastoid';
                            %                     ElectrodeManufacturer='PMT';
                            iEEGGround='chest';
                            iEEGPlacementScheme='clinical indication';
                            iEEGReference=ChannelInformation{5,16};
                            ElectrodeManufacturer=ChannelInformation{1,6};
                            GroupList='multiple sEEG depths, see electrode locations in electrodes.tsv file';
                            MRIDataAlignedTo='defacedMRI';
                            
                            BIPOL=0;
                            % you can also have acq- and proc-, but these are optional
                            %% Preimplant scan folder and copy file. These data should be defaced.
                            preimpfolder = fullfile(root_dir, project_label, ...
                                ['sub-' ieeg_sub], ['ses-' 'preimp']);
                            preimpMRI = fullfile(root_dir, project_label, ...
                                ['sub-' ieeg_sub], ['ses-' 'preimp'], ['anat']);
                            %                            if czv==1
                            acq_id = 'T1w';
                            ieeg_sesanat='preimp';
                            
                            %A json file is needed to describe the scan and
                            %the information
                            createBIDS_T1w_json_short
                            
                            preimMRIName = fullfile(['sub-' ieeg_sub ...
                                '_ses-' ieeg_sesanat ...
                                '_acq-' acq_id ...
                                '_run-' '01' '_T1w.nii']);
                            
                            mkdir(preimpMRI)
                            copyfile([Maindir,PatientName,'\',MRIDir(1).name],[preimpMRI,'\',preimMRIName])
                            
                            %% Post implant scan which can be CT or MRI. Requires the same steps as the preimplant scan. These data should be defaced.
                            postimpscan = fullfile(root_dir, project_label, ...
                                ['sub-' ieeg_sub], ['ses-' 'postimp'], ['anat']);
                            if DI<72
                                %                                 postimMRIName=fullfile(['sub-' ieeg_sub ...
                                %                                     '_ses-' ieeg_ses ...
                                %                                     '_CT.nii']);
                                acq_id = 'CT';
                                ieeg_sesanat='postimp';
                                createBIDS_T1w_json_short
                                postimMRIName = fullfile(['sub-' ieeg_sub ...
                                    '_ses-' ieeg_sesanat ...
                                    '_acq-' acq_id ...
                                    '_run-' '01' '_T1w.nii']);
                            else
                                %                                 postimMRIName=fullfile(['sub-' ieeg_sub ...
                                %                                     '_ses-' ieeg_ses ...
                                %                                     '_T1w.nii']);
                                acq_id = 'T1w';
                                createBIDS_T1w_json_short
                                postimMRIName = fullfile(['sub-' ieeg_sub ...
                                    '_ses-' ieeg_sesanat ...
                                    '_acq-' acq_id ...
                                    '_run-' '01' '_T1w.nii']);
                            end
                            mkdir(postimpscan)
                            copyfile([Maindir,PatientName,'\',PostopMRIDir(1).name],[postimpscan,'\',postimMRIName])
                            
                            
                            %% This is where the ieeg data per task/session data along with the channel information is organized and saved:
                            MainDirBIDS = fullfile(root_dir, project_label, ...
                                ['sub-' ieeg_sub], ['ses-' ieeg_ses]);
                            
                            postimpFolder = fullfile(root_dir, project_label, ...
                                ['sub-' ieeg_sub], ['ses-' 'postimp']);
                            mkdir([postimpFolder,'\ieeg'])
                            
                            createBIDS_channels_tsvPerBank
                            
                            createBIDS_electrodes_tsvStandardization
                            
                            createBIDS_coordsystem_json2
                            
                            createBIDS_ieeg_jsonFileOnly
                            
                            createBIDS_events_tsv_json_full
                            
                            createBIDS_data_WriteBrainVisionWithFieldtrip
                            
                            SummaryArray(DI,1)=Patients(DI);
                            SummaryArray(DI,2)=Relabel(DI);
                            SummaryArray(DI,3)=Identifiers(DI);
                            if NT==1
                                SummaryArray(DI,czv+3)={'Single Pulse, Multiple stimulation sites'};
                            elseif NT==2
                                SummaryArray(DI,czv+3)={'Single Pulse, Multiple stimulation current amplitudes'};
                            elseif NT==3
                                SummaryArray(DI,czv+3)={'Trains, Multiple frequencies and Current Amplitudes, Single pulse'};
                            elseif NT==4
                                SummaryArray(DI,czv+3)={'Single Pulse, Multiple stimulation sites'};
                            elseif NT==5
                                SummaryArray(DI,czv+3)={'Single Pulse, Multiple stimulation sites'};
                            end
                            
                            czv=czv+1
                        end
                    end
                end
            end
        end
        MatchedIDs(DI,1)=Patients(DI);
        MatchedIDs(DI,2)=Identifiers(DI);
        MatchedIDs(DI,3)=Relabel(DI);
    end
end



% writecell(SummaryArray,'Y:\StimDataBackup\Data_Stimulation\DeIdentifiedDataSet\test.xls','Sheet',1)

