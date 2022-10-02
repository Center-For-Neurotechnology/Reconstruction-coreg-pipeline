% clear
% addpath(genpath('X:\Projects\Lab_Materials\Analysis_Tools_and_Software\fieldtrip-20220202\'))
addpath(genpath('Y:\ReconPipelinePaper\Code\'))

SubcorticalLabels={'Brain-Stem';'Left-Accumbens-area';'Left-Amygdala';...
    'Left-Caudate';'Left-Cerebellum-Cortex';'Left-Hippocampus';...
    'Left-Pallidum';'Left-Putamen';'Left-Thalamus-Proper';'Left-Thalamus';....
    'Right-Accumbens-area';'Right-Amygdala';'Right-Caudate';....
    'Right-Cerebellum-Cortex';'Right-Hippocampus';'Right-Pallidum';....
    'Right-Putamen';'Right-Thalamus-Proper';'Right-Thalamus'};

Patients{1}='sub-0t3i';

MissingFile=[];MissingFileLabel={};MissingFileSubcort=[];MissingFileLabelSubcort={};
for DI=1:size(Patients,1)
    %     if isempty(dir(['Y:\StimDataBackup\ReconLocs\',Patients{DI},'\',Patients{DI},'MappedRegionsVertices*']))==1
    Patients{DI}
    PtName=Patients{DI};
    MRIDirectory=['Y:\ReconPipelinePaper\Data\derivatives\freesurfer\',PtName,'\surf\'];
    RASExcelDirectory=['Y:\ReconPipelinePaper\Data\',PtName,'\ses-postimp\ieeg\'];
    ELAParcelDirectory=['Y:\ReconPipelinePaper\Data\derivatives\MMVT\',PtName,'\electrodes\'];
    VolumeDir=['Y:\ReconPipelinePaper\Data\derivatives\ParcellatedVolumes\',PtName,'\'];
    
    if isempty(dir([RASExcelDirectory,'\',PtName,'_ElectrodeInfo_RAS.xlsx*']))==0 &&...
            isempty(dir([MRIDirectory,'*rh.pial']))==0
        %% Creates vertices in 3D that can be used to measure the nearest pial surface, white matter boundary, or subcortical region
        
        [verticesrh, facesrh] = freesurfer_read_surf([MRIDirectory,'rh.pial']);
        [verticesWhiterh, facesWhiterh] = freesurfer_read_surf([MRIDirectory,'rh.white']);
        
        [verticeslh, faceslh] = freesurfer_read_surf([MRIDirectory,'lh.pial']);
        [verticesWhitelh, facesWhitelh] = freesurfer_read_surf([MRIDirectory,'lh.white']);
        
        SurfOuter=[verticesrh;verticeslh];
        SurfInner=[verticesWhiterh;verticesWhitelh];
        
        SCDirbase=dir([VolumeDir,'\aparc_',SubcorticalLabels{1},'*.stl']);
        ValVertSubcort=[];SCLabel={};ScVertFace=[];
        if isempty(SCDirbase)==0
            for SC=1:length(SubcorticalLabels)
                SCDir=dir([VolumeDir,'\aparc_',SubcorticalLabels{SC},'*.stl']);
                if isempty(SCDir)==0
                    FileLoad=[VolumeDir,SCDir(1).name];
                    gm = stlread(FileLoad);
                    shp=alphaShape(gm.Points);
                    ValVertSubcort=[ValVertSubcort;repmat(SC,size(gm.Points,1),1) gm.Points];
                    ScVertFace(SC).faces=shp;
                    ScVertFace(SC).vertices=gm.Points;
                    %             ValVertSubcort=[ValVertSubcort;SC];
                    SC
                    SCLabel{SC}=SCDir(1).name;
                else
                    ScVertFace(SC).faces=[];
                    ScVertFace(SC).vertices=[];
                    SCLabel{SC}='';
                end
            end
        end
        
        %% Reads the RAS spreadsheet
        RASTable=readtable([RASExcelDirectory,'\',PtName,'_ElectrodeInfo_RAS.xlsx']);
        RAS_coords=table2array(RASTable(:,2:4));
        RAS_labels=table2array(RASTable(:,1));
        RAS=table2array(RASTable(:,2:4));
        RASlabel=table2array(RASTable(:,1));
        ElectrodeType=table2array(RASTable(:,5));
        Manufacturer=table2array(RASTable(:,6));
        SurfaceAreasize=table2array(RASTable(:,7));
        group=table2array(RASTable(:,8));
        hemisphere=table2array(RASTable(:,9));
        R=RAS(:,1);
        A=RAS(:,2);
        S=RAS(:,3);
        
        Chankeep=find(isnan(RAS_coords(:,1))==0);
        RAS_coords=RAS_coords(Chankeep,:);
        RAS_labels=RAS_labels(Chankeep);
        if size(RAS_coords,2)>3
            RAS_coords=RAS_coords(:,1:3);
        end
        
        %%                 Getting a shape for the outer pial surface and the white
        %                 matter to identify contacts in different volumes.
        shpOuter = alphaShape(SurfOuter(:,1),SurfOuter(:,2),SurfOuter(:,3));
        shpInner = alphaShape(SurfInner(:,1),SurfInner(:,2),SurfInner(:,3));
        
        shpSC = alphaShape(ValVertSubcort(:,2),ValVertSubcort(:,3),ValVertSubcort(:,4),10);
        shpInnerConv = alphaShape(SurfInner(:,1),SurfInner(:,2),SurfInner(:,3),Inf);
        
        %% Reads per-contact and identifies if the contact is in grey matter, white matter, or outside brain or subcortex
        %             As well as measures the distance to the closest grey matter,
        %             white matter boundary or subcortical region
        
        
        RASDistPerLoc=[];RASClosestPerLoc=[];
        
        NearestWhiteMatterTract=cell(length(RAS_coords),1);
        NearestWhiteMatterTractIndex=zeros(length(RAS_coords),1);
        
        SingleCoordLoc=[];
        DistanceGreyWhite=NaN*ones(length(RAS_coords),1);
        DistancePial=NaN*ones(length(RAS_coords),1);
        NumberOfContactsNearby20mm=NaN*ones(length(RAS_coords),1);
        for Electrode=1:length(RAS_coords)
            
            Coords=[RAS_coords(Electrode,:)];
            
            %% Identifying if the single contact is in grey matter, white matter, or subcortical.
            inWM = inShape(shpInner,Coords(:,1),Coords(:,2),Coords(:,3));
            inGM = inShape(shpOuter,Coords(:,1),Coords(:,2),Coords(:,3));
            inSC = inShape(shpSC,Coords(:,1),Coords(:,2),Coords(:,3));
            inWMConv = inShape(shpInnerConv,Coords(:,1),Coords(:,2),Coords(:,3));
            
            if inWM==1 && inSC==1 && inGM==1 && inWMConv==1 %Subcortical
                SingleCoordLoc=[SingleCoordLoc;Coords 3 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 3 Electrode Electrode inWM inGM inSC];
            elseif inWM==1 && inSC==0 && inGM==1 && inWMConv==1 %grey-white matter overlap
                SingleCoordLoc=[SingleCoordLoc;Coords 5 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 5 Electrode Electrode inWM inGM inSC];
            elseif inWM==0 && inSC==0 && inGM==1 && inWMConv==0 %grey matter
                SingleCoordLoc=[SingleCoordLoc;Coords 2 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 2 Electrode Electrode inWM inGM inSC];
            elseif inWM==0 && inSC==1 && inGM==1 && inWMConv==1 %grey matter and subcortical
                SingleCoordLoc=[SingleCoordLoc;Coords 6 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 6 Electrode Electrode inWM inGM inSC];
            elseif inWM==1 && inSC==0 && inGM==0 && inWMConv==1 %white matter
                SingleCoordLoc=[SingleCoordLoc;Coords 4 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 4 Electrode Electrode inWM inGM inSC];
            elseif inWM==0 && inSC==0 && inGM==0 && inWMConv==1 %white matter
                SingleCoordLoc=[SingleCoordLoc;Coords 4 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 4 Electrode Electrode inWM inGM inSC];
            elseif inWM==0 && inSC==1 && inGM==0 && inWMConv==1 %subcortical
                SingleCoordLoc=[SingleCoordLoc;Coords 3 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 3 Electrode Electrode inWM inGM inSC];
            elseif inWM==0 && inSC==0 && inGM==0 && inWMConv==0 %outside brain
                SingleCoordLoc=[SingleCoordLoc;Coords 1 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 1 Electrode Electrode inWM inGM inSC];
            elseif inWM==0 && inSC==0 && inGM==1 && inWMConv==1 %grey matter
                SingleCoordLoc=[SingleCoordLoc;Coords 2 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 2 Electrode Electrode inWM inGM inSC];
            else %other
                SingleCoordLoc=[SingleCoordLoc;Coords 7 Electrode Electrode inWM inGM inSC];
                CoorPerLoc(Electrode,:)=[Coords 7 Electrode Electrode inWM inGM inSC];
            end
            
            LocMap={'outside brain','grey matter','subcortical',...
                'white matter','grey-white matter overlap','grey matter and subcortical','other'};
            
            %% distances to the outer pial surface, the nearest grey-white boundary, and the nearest subcortical region
            %compute Euclidean distances to the outer pial surface:
            distancesGM = sqrt(sum(bsxfun(@minus, SurfOuter, Coords).^2,2));
            closestGM = SurfOuter(distancesGM==min(distancesGM),:);
            distclosestGM = distancesGM(distancesGM==min(distancesGM),:);
            %                 v1=CoordsBipol-repmat(CoordsBipol(1,:),2,1);
            %                 v2=closest-repmat(CoordsBipol(1,:),1,1);
            %                 [magn, angle] = ab2v(v1(2,:), v2(1,:));
            
            %compute Euclidean distances  to the white matter boundary:
            distancesWM = sqrt(sum(bsxfun(@minus, SurfInner, Coords).^2,2));
            closestWM = SurfInner(distancesWM==min(distancesWM),:);
            distclosestWM = distancesWM(distancesWM==min(distancesWM),:);
            %                 v1=Coords-repmat(Coords(1,:),2,1);
            %                 v2=closestWM-repmat(Coords(1,:),1,1);
            %                 [magnWM, angleWM] = ab2v(v1(2,:), v2(1,:));
            
            DistanceGreyWhite(Electrode)=distclosestGM-distclosestWM;
            DistancePial(Electrode)=distclosestGM;
            
            distancesWMGM=[];closestWMGM=[];
            if isempty(closestWM)==0
                distancesWMGM = sqrt(sum(bsxfun(@minus, SurfOuter, closestWM).^2,2));
                closestWMGM = SurfOuter(distancesWMGM==min(distancesWMGM),:);
                %                     v1=Coords-repmat(Coords(1,:),2,1);
                %                     v2=closestWMGM-repmat(Coords(1,:),1,1);
                %                     [magnWMGM, angleWMGM] = ab2v(v1(2,:), v2(1,:));
            end
            distancesSC=[];closestSC=[];
            if isempty(ValVertSubcort)==0
                %compute Euclidean distances to subcortical structures:
                distancesSC = sqrt(sum(bsxfun(@minus, ValVertSubcort(:,2:end), Coords).^2,2));
                %find the smallest distance and use that as an index into B:
                closestSC = ValVertSubcort(distancesSC==min(distancesSC),2:4);
                %                     v1=Coords-repmat(Coords(1,:),2,1);
                %                     v2=closestSC-repmat(Coords(1,:),1,1);
                %                     [magnSC, angleSC] = ab2v(v1(2,:), v2(1,:));
            else
                MissingFileSubcort(DI,:)=[DI isempty(dir([RASExcelDirectory,'\',PtName,'_ElectrodeInfo_RAS.xlsx*'])) isempty(dir([MRIDirectory,'*rh.pial'])) isempty(dir([aparc_',SubcorticalLabels{1},'*.stl']))];
                MissingFileLabelSubcort{DI}=PtName;
            end
            
            %% Nearest neighbor electrodes (within 20 mm)
            distancesOtherElec = sqrt(sum(bsxfun(@minus, RAS_coords, Coords).^2,2));
            closestElecLoc = RAS_coords(distancesOtherElec<20,:);
            closestElecLabel = RAS_labels(distancesOtherElec<20);
            
            NumberOfContactsNearby20mm(Electrode)=length(closestElecLabel);
            %% Saved structure with detailed information
            
            RASClosestPerLoc(Electrode).closest=closestGM;
            RASClosestPerLoc(Electrode).closestWM=closestWM;
            RASClosestPerLoc(Electrode).closestSC=closestSC;
            RASClosestPerLoc(Electrode).Coords=Coords;
            RASClosestPerLoc(Electrode).distancesOtherElec=distancesOtherElec;
            RASClosestPerLoc(Electrode).closestElecLoc=closestElecLoc;
            RASClosestPerLoc(Electrode).closestElecLabel=closestElecLabel;
            RASClosestPerLoc(Electrode).WhiteToGreyCoord=closestWMGM;
            
            NearestWhiteMatterTract{Electrode}='unknown';
            NearestWhiteMatterTractIndex(Electrode)=0;
        end
        
        NearGreyWhiteBoundary=(CoorPerLoc(:,4)>4 & CoorPerLoc(:,4)<7);
        InsideBrain=CoorPerLoc(:,4)~=1;
        VolumeLabelIndex=CoorPerLoc(:,4);
        VolumeLabel=LocMap(VolumeLabelIndex)';
        
        %% Per electrode depth information
        
        CL=char(RAS_labels);
        LabelOnly={};
        ChannelNumber=[];
        for lc=1:length(CL)
            Nm=CL(lc,:);
            NumberInd=[];
            for ms=1:length(Nm)
                if isempty(str2num(Nm(ms)))==0
                    NumberInd=[NumberInd ms];
                end
            end
            if isempty(NumberInd)==1
                LabelOnly{lc}=Nm;
                ChannelNumber(lc)=1;
            else
                if isempty(findstr('MG72',PtName))==0 && isempty(findstr('L',Nm))==0
                    LabelOnly{lc}=Nm(1:NumberInd(2)-1);
                    ChannelNumber(lc)=str2num(Nm(NumberInd(2):end));
                else
                    LabelOnly{lc}=Nm(1:NumberInd(1)-1);
                    ChannelNumber(lc)=str2num(Nm(NumberInd(1):end));
                end
            end
        end
        
        %% ELA localization
        MappingELAParcellationsScript
        
        ProbabilityMapping(isnan(ProbabilityMapping(:,8)),8)=length(TargetLabelsParc);
        ParcellationELALabel=TargetLabelsParc(ProbabilityMapping(:,8));
        ParcellationELALabelIndex=ProbabilityMapping(:,8);
        ProbabilityELAGrey=ProbabilityMapping(:,2);
        ProbabilityELAWhite=ProbabilityMapping(:,3);
        
        %% Measures along a single depth or grid (grouped by per depth or per grid)
        [INElec,ElecList]=grp2idx(LabelOnly);
        SingleCoordLoc=[];LeadDist=[];
        InterContactDistBipol=[];
        for Devices=1:length(ElecList)
            
            IndCheck=find(INElec==Devices);
            CoordsLength=([RAS_coords(IndCheck,:)]);
            PM=ProbabilityMapping(IndCheck,:);
            distancesperlead=[];
            if size(CoordsLength,1)>1
                for coor2=1:size(CoordsLength,1)-1
                    distancesperlead(coor2) = sqrt(sum(bsxfun(@minus, CoordsLength(coor2,:), CoordsLength(coor2+1,:)).^2,2));
                end
            end
            
            InterContactDistBipol=[InterContactDistBipol; [NaN; distancesperlead']];
            
            DistLead=sqrt((CoordsLength(end,1)-CoordsLength(1,1))^2+...
                (CoordsLength(end,2)-CoordsLength(1,2))^2+...
                (CoordsLength(end,3)-CoordsLength(1,3))^2);
            
            LeadDist(Devices).distancesperlead=distancesperlead;
            LeadDist(Devices).DistLead=DistLead;
            LeadDist(Devices).CoordsLen=size(CoordsLength,1);
            LeadDist(Devices).Coords=CoordsLength;
        end
        
        %% Volumetric parcellation
        
        Volumes=dir([VolumeDir,'*.stl']);
        
        ShapeOverall=[];VolHeaderLabel={};LabelVolume=[];
        LabelWhiteMatter=[];LabelCortex=[];
        for FiCol=1:length(Volumes)
            FileLoad=[VolumeDir,Volumes(FiCol).name];
            
            %read the stl as a triangulation
            gm = stlread(FileLoad);
            
            shp=alphaShape(gm.Points);
            inshp= inShape(shp,RAS(:,1),RAS(:,2),RAS(:,3));
            ShapeOverall(:,FiCol)=inshp;
            sfwe=1;
            LabelAtlas=Volumes(FiCol).name(1:end-4);
            LabelAtlas=(LabelAtlas(7:end));
            VolHeaderLabel{FiCol}=LabelAtlas;
            LabelVolume(FiCol).triangulationShape=gm;
            LabelVolume(FiCol).Vertices=gm.Points;
            LabelVolume(FiCol).AlphaShape=shp;
            LabelVolume(FiCol).LabelAtlas=LabelAtlas;
            if contains(LabelAtlas,'Cerebral-White')==1
                LabelWhiteMatter=[LabelWhiteMatter;FiCol];
            end
            if contains(LabelAtlas,'Cerebral-Cortex')==1 || contains(LabelAtlas,'unknown')==1
                LabelCortex=[LabelCortex;FiCol];
            end
            
            FiCol
        end
        VolHeaderLabel{length(Volumes)+1}='unknown';
        
        ShapeOverallLabels=ShapeOverall;
        ShapeOverallwmcX=ShapeOverall;
        
        ShapeOverallLabels(:,[LabelCortex;LabelWhiteMatter])=0;
        
        [mx,ind]= max(ShapeOverallLabels');
        [mxwmcX,indwmcX]= max(ShapeOverallwmcX');
        
        indwmcX(mxwmcX==0)=length(Volumes)+1;
        ind(mx==0)=indwmcX(mx==0);
        
        ParcellationVolLabel=VolHeaderLabel(ind)';
        ParcellationVolLabelIndex=ind';
        
        T = table(RAS_labels,R,A,S,ElectrodeType,Manufacturer,...
            SurfaceAreasize,group,hemisphere,...
            ParcellationVolLabel,ParcellationVolLabelIndex,...
            ParcellationELALabel,ParcellationELALabelIndex,...
            ProbabilityELAGrey,ProbabilityELAWhite,...
            InsideBrain,DistanceGreyWhite,DistancePial,VolumeLabel,VolumeLabelIndex,...
            InterContactDistBipol,NearGreyWhiteBoundary,...
            NumberOfContactsNearby20mm,...
            NearestWhiteMatterTract,NearestWhiteMatterTractIndex);
        
        writetable(T,['Y:\ReconPipelinePaper\Data\derivatives\',PtName,'_RAS_Electrode_Parc.csv'])
        writetable(T,['Y:\ReconPipelinePaper\Data\derivatives\',PtName,'_RAS_Electrode_Parc.xlsx'])
        
        BrainLocationLabelsELA=BrainLocationLabels;
        save(['Y:\ReconPipelinePaper\Data\derivatives\',PtName,'_RAS_Electrode_Parc'],...
            'T','LabelVolume','VolHeaderLabel','RASClosestPerLoc',...
            'RAS_coords','RAS_labels',...
            'ValVertSubcort','SCLabel','facesWhiterh','facesWhitelh','verticesWhiterh',...
            'verticesWhitelh','SurfInner','SurfOuter','verticesrh',...
            'verticeslh','faceslh','facesrh','CoorPerLoc','ScVertFace','LeadDist',...
            'AllProbabilitiesPerSite','BrainLocationLabelsELA','TargetLabelsParc')
        
        
    else
        MissingFile(DI,:)=[DI isempty(dir([RASExcelDirectory,'\',PtName,'_ElectrodeInfo_RAS.xlsx*'])) isempty(dir([MRIDirectory,'*rh.pial'])) isempty(dir([VolumeDir,'*']))];
        MissingFileLabel{DI}=PtName;
    end
    %     end
end
%%
BrainCol=[.95 .95 .95];
COLVol=colormap(hsv(length(LocMap)));
COLVolReg=colormap(hsv(length(VolHeaderLabel)));
COLELAReg=colormap(hsv(length(TargetLabelsParc)));
colormap(hsv(255))
clf
close all
figure(1)
subplot(2,2,1)
patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
hold on
patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
axis off
set(gcf,'color','w')
scatter3(RAS_coords(:,1), RAS_coords(:,2), RAS_coords(:,3),8,[0 0 0],'filled')
scatter3(RAS_coords(InsideBrain==0,1), RAS_coords(InsideBrain==0,2), RAS_coords(InsideBrain==0,3),32,'r','filled')
text(RAS_coords(isnan(InterContactDistBipol),1), RAS_coords(isnan(InterContactDistBipol),2), RAS_coords(isnan(InterContactDistBipol),3),RAS_labels(isnan(InterContactDistBipol)))
title('Electrodes in brain (electrodes outside brain in red)')
daspect([1 1 1])

% figure(2)
subplot(2,2,2)
patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
hold on
patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
axis off
set(gcf,'color','w')
% scatter3(RAS_coords(:,1), RAS_coords(:,2), RAS_coords(:,3),8,[0 0 0],'filled')
for kw=1:length(LocMap)
    plot3(RAS_coords(VolumeLabelIndex==kw,1), RAS_coords(VolumeLabelIndex==kw,2),...
        RAS_coords(VolumeLabelIndex==kw,3),'color',COLVol(kw,:),'markersize',12,...
        'linestyle','none','marker','.')
end
LOCS={{'left H'},{'right H'},LocMap}
legend(cat(1,{'left H'},{'right H'},LocMap{:}))
text(RAS_coords(isnan(InterContactDistBipol),1), RAS_coords(isnan(InterContactDistBipol),2), RAS_coords(isnan(InterContactDistBipol),3),RAS_labels(isnan(InterContactDistBipol)))
title('Electrodes in white matter, grey matter, and other categories')
daspect([1 1 1])

% figure(3)
subplot(2,2,3)
patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
hold on
patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
axis off
set(gcf,'color','w')
% scatter3(RAS_coords(:,1), RAS_coords(:,2), RAS_coords(:,3),8,[0 0 0],'filled')
for kw=1:length(VolHeaderLabel)
    plot3(RAS_coords(ParcellationVolLabelIndex==kw,1), RAS_coords(ParcellationVolLabelIndex==kw,2),...
        RAS_coords(ParcellationVolLabelIndex==kw,3),'color',COLVolReg(kw,:),'markersize',12,...
        'linestyle','none','marker','.')
end
LOCS={{'left H'},{'right H'},VolHeaderLabel}
% legend(cat(1,{'left H'},{'right H'},VolHeaderLabel{:}))
title('Electrode locations, Volume mapped')
daspect([1 1 1])

% figure(4)
subplot(2,2,4)
patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
hold on
patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[BrainCol],'FaceAlpha',.1,'EdgeColor','none')
axis off
set(gcf,'color','w')
for kw=1:length(TargetLabelsParc)
    plot3(RAS_coords(ParcellationELALabelIndex==kw,1), RAS_coords(ParcellationELALabelIndex==kw,2),...
        RAS_coords(ParcellationELALabelIndex==kw,3),'color',COLELAReg(kw,:),'markersize',12,...
        'linestyle','none','marker','.')
end
LOCS={{'left H'},{'right H'},TargetLabelsParc}
% legend(cat(1,{'left H'},{'right H'},TargetLabelsParc{:}))

title('Electrode locations, ELA')
daspect([1 1 1])


