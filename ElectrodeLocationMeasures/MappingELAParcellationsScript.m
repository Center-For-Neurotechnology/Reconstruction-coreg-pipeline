clf
COdi=1;

%add precentral and postcentral
co=1;

ParcellationValues=NaN*ones(size(RAS_coords,1),8);

% This is the targeted labels I start with to be able to classify the
% channels.
TargetLabelsParc={'rostralmiddlefrontal'
    'superiorfrontal'
    'caudalmiddlefrontal'
    'medialorbitofrontal'
    'lateralorbitofrontal'
    'parstriangularis'
    'accumbens'
    'parsorbitalis'
    'parsopercularis'
    'superiortemporal'
    'middletemporal'
    'inferiortemporal'
    'transversetemporal'
    'rostralanteriorcingulate'
    'caudalanteriorcingulate'
    'posteriorcingulate'
    'amygdala'
    'hippocamp'
    'parahippocampal'
    'insula'
    'caudate'
    'precentral'
    'postcentral'
    'inferiorparietal'
    'superiorparietal'
    'fusiform'
    'putamen'
    'postcentral'
    'thalamus'
    'lingual'
    'pericalcarine'
    'occipital'
    'isthmuscingulate'
    'supramarginal'
    'precuneus'
    'cuneus'
    'cerebellum'
    'ventraldc'
    'pallidum'
    'paracentral'
    'entorhinal'
    'CC_Anterior'
    'hypointensities'
    'NotLabelled'
    'unknown'};

PatientName=PtName;
ParcelFile=[PtName,'_aparc.DKTatlas40_electrodes_cigar_r_3_l_4.csv'];

ParcelELA1=readtable([ELAParcelDirectory,'\',ParcelFile]);
    TargetLabelsParc={};
    TargetLabelsParc=ParcelELA1.Properties.VariableNames(2:end)';
    
TimeMin=-1;
valsig=10;
subtTime=1;

TASK1='Amplitude';
Rec=[];

MissingFile=[];MissingFileLabel={};MissingFileSubcort=[];MissingFileLabelSubcort={};NotMissingFileSubcort=[];NotMissingFileLabelSubcort={};

ParcDir=dir([ELAParcelDirectory,'\*',ParcelFile]);

if isempty(ParcDir)==0
    %Parcellation file loading and organization
    ParcelELA=readtable([ELAParcelDirectory,'\',ParcelFile]);
    ElectLabelRAS=table2array(ParcelELA(:,1));
    ElecLocsParc={};
    ElecLocsParc=ElectLabelRAS;
    ParcellationMatrix=table2array(ParcelELA(:,2:end));
    ParcellationMatrix2=ParcellationMatrix(isnan(ParcellationMatrix(:,1))==0,:);
    BrainLocationLabels={};
    BrainLocationLabels=ParcelELA.Properties.VariableNames(2:end);
    BrainLocationProbabilities=[];
    BrainLocationProbabilities=ParcellationMatrix2;
    
    WM=[];elec=[]; approx=[];WMR=[];WML=[];
    for MIN=1:length(BrainLocationLabels)
        if contains(BrainLocationLabels{MIN},'Cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'Right','IgnoreCase',true)==1
            WMR=[MIN];
        elseif contains(BrainLocationLabels{MIN},'Cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'Left','IgnoreCase',true)==1
            WML=[MIN];
        end
        if contains(BrainLocationLabels{MIN},'Cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'rh','IgnoreCase',true)==1
            WMR=[MIN];
        elseif contains(BrainLocationLabels{MIN},'Cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'lh','IgnoreCase',true)==1
            WML=[MIN];
        end
        if contains(BrainLocationLabels{MIN},'cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'Right','IgnoreCase',true)==1
            WMR=[MIN];
        elseif contains(BrainLocationLabels{MIN},'cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'Left','IgnoreCase',true)==1
            WML=[MIN];
        end
        if contains(BrainLocationLabels{MIN},'cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'rh','IgnoreCase',true)==1
            WMR=[MIN];
        elseif contains(BrainLocationLabels{MIN},'cerebral','IgnoreCase',true)==1 && contains(BrainLocationLabels{MIN},'lh','IgnoreCase',true)==1
            WML=[MIN];
        end
        if contains(BrainLocationLabels{MIN},'elc_length','IgnoreCase',true)==1
            elec=[MIN];
        end
        if contains(BrainLocationLabels{MIN},'approx','IgnoreCase',true)==1
            approx=[MIN];
        end
    end
    
    BrainLocationProbabilities(:,[WMR WML elec approx])=0;

    CL=char(RAS_labels);
    [IM,Label]=grp2idx(CL(:,1:3));
    [MI,NonLRelecLabels]=grp2idx(CL(:,2:3));
    TargetChannels=[];
    COL=colormap(cool(max(MI)));
    
    %This steps through all the RAS channels, gets the
    %locations for a bipolar pair, and matches it to the
    %label probabilities and white matter probabilities.
    ProbAll=[];ProbWM=[]; RASLocLabel=[];LocLabel=[];cou=1;ProbEle=[];Gonogo=0;
    ProbabilityMapping=[];
    AllProbabilitiesPerSite=[];
    ProbabilityMapping=NaN*ones(size(RAS_coords,1),8);
    for Ele=1:length(ElecLocsParc)
        MainElec=ElecLocsParc{Ele};
        ElecLocTarg= MainElec;
        if isempty(strmatch(MainElec(1),'R'))==0
            WM=WMR;
            WhiteMatterProb= ParcellationMatrix2(Ele,WM)+10^-6;
        elseif isempty(strmatch(MainElec(1),'L'))==0
            WM=WML;
            WhiteMatterProb= ParcellationMatrix2(Ele,WM)+10^-6;
        end
        if isempty(WM)==1
            WM=WML;
            WhiteMatterProb= max(ParcellationMatrix2(Ele,[WML WMR]))+10^-6;
        end
        
        [MX,IX]= max(BrainLocationProbabilities(Ele,:));
        if MX>0
            BrainReg=BrainLocationLabels{IX};
            Elecprob=MX;
        else
            BrainReg='unknown';
            Elecprob=MX;
        end
        if length(find(BrainLocationProbabilities(Ele,:)>0))>1
            AllProbabilitiesPerSite(Ele).Labels=BrainLocationLabels(BrainLocationProbabilities(Ele,:)>0);
            AllProbabilitiesPerSite(Ele).Indices=find(BrainLocationProbabilities(Ele,:)>0);
            AllProbabilitiesPerSite(Ele).Elecprobability=BrainLocationProbabilities(Ele,BrainLocationProbabilities(Ele,:)>0);
        end
        
        ProbAll1=log(Elecprob/WhiteMatterProb);
        ProbWM1=(WhiteMatterProb);
        ProbEle1=(Elecprob);
        cou=1;
        %%
        Gonogo=0;LabelValue=NaN; ProbEle=NaN;  ProbWM=NaN; ProbAll=NaN; RASLocLabel=[NaN NaN NaN];
        for MN=1:length(CL)
            ELP=ElecLocTarg;
            ELP(findstr('-',ELP))=' ';
            Chan1=[ELP,' '];
            ChanRASList1=CL(MN,:);
            if isempty(findstr('0',ChanRASList1(4)))==0
                ChanRASList1(findstr('0',ChanRASList1))=[];
            end
            if isempty(findstr(ChanRASList1,Chan1))==0
                LocLabel=ChanRASList1;
                RASLocLabel=nanmean(RAS_coords(MN,:),1);
                ProbAll=log(Elecprob/WhiteMatterProb);
                ProbWM=(WhiteMatterProb);
                ProbEle=(Elecprob);
                
                Gonogo=0;LabelValue=NaN;
                for LH=1:length(TargetLabelsParc)
                    if isempty(findstr(TargetLabelsParc{LH},BrainReg))==0
                        Gonogo=1;
                        LH
                        
                        LabelValue=LH;
                    end
                end
            end
        end
        ProbabilityMapping(Ele,:)=[ Gonogo  ProbEle  ProbWM ProbAll RASLocLabel LabelValue];
    end
    
    ProbabilityMappingHeader={'Gonogo(found a matching parcellation label)', 'Probability of belonging to that label','Probability of belonging in white matter',...
        'log(label probability/white matter probability)','RAS coordinate 1','RAS coordinate 2',...
        'RAS coordinate 3','which label in the TargetLabels list of regions'};
    
end

