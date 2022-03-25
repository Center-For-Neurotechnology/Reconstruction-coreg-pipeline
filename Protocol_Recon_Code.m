%% Coregistration
% This Code is the Forerunner of All Good Reconstruction Code
% In Development: Jan 2022
% Note To Self: Run All Previous Sections Before Running Last Section


PatientName='id##';
wd = '/home/bourbon-the-huckster/Documents/Recons';

DirVal = fullfile(wd , PatientName, [PatientName '_SurferOutput'], 'surf'); %the surf folder with the freesurfer output
ImageDirectory = fullfile(wd , PatientName, 'Images');         %print files
RASExcelDirectory = fullfile(wd , PatientName);
MRIDirectory = fullfile(wd , PatientName, [PatientName '_SurferOutput'], 'mri'); %mri.nii...NOT nii.gz

[num,txt,all]=xlsread( fullfile( RASExcelDirectory, [PatientName '_RAS.xlsx'] ) );
RAS_coords=[];
RAS_coords=num;
RAS_labels={};
RAS_labels=txt(2:end,1);
% RAS_coords_fixed=[];
% RAS_no_depths=[];

% Check for Images Directory in the Working Directory, make it
% if it's not there
if ~isfolder(ImageDirectory)
    
    mkdir(ImageDirectory)
    disp('Made the Images folder in your Patient''s Working Directory!')
    
end


csvwrite([RASExcelDirectory, PatientName,'_RAS.csv'],RAS_coords)
display("Made it Through. Check the Workspace, that the files might match your expectations.");

%% Getting the freesurfer pial and seeing where things are located, with the white matter in blue and pial surface.

% It's also handy for the blender-friendly stl file saving.
[verticesrh, facesrh] = freesurfer_read_surf( fullfile(DirVal,'rh.pial') );
stlwrite( fullfile(DirVal,'right.stl'), facesrh, verticesrh)

[verticeslh, faceslh] = freesurfer_read_surf( fullfile(DirVal,'lh.pial') );
stlwrite( fullfile(DirVal,'left.stl'), faceslh, verticeslh)

[verticesrhwhite, facesrhwhite] = freesurfer_read_surf( fullfile(DirVal,'rh.white') );
stlwrite( fullfile(DirVal,'rightWhite.stl'), facesrhwhite, verticesrhwhite)

[verticeslhwhite, faceslhwhite] = freesurfer_read_surf( fullfile(DirVal,'lh.white') );
stlwrite([DirVal,'leftWhite.stl'], faceslhwhite, verticeslhwhite)

%Surfer Output Tester
patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[.7 .7 .7],'facealpha',.9,'edgecolor','none')
hold on
patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[.7 .7 .7],'facealpha',.9,'edgecolor','none')
camlight('headlight');

h=gcf;
[cortexr.vert, cortexr.tri] = read_surf( fullfile(DirVal, 'lh.pial') );
[cortexl.vert, cortexl.tri] = read_surf( fullfile(DirVal, 'rh.pial') );


tripatch(cortexr,h,[.5 .5 .5],'facealpha',1);
hold on
tripatch(cortexl,h,[.5 .5 .5],'facealpha',1);
camlight('headlight');
shading interp; lighting gouraud; material dull;


%% Plotting RAS electrodes on the Surfer Output
clf

% Create figure handle in fig
fig=gcf;

%Returns a matrix of vertex coordinates and face lists for the given .pial file
[cortexr.vert, cortexr.tri] = read_surf( fullfile(DirVal, 'rh.pial') );
[cortexl.vert, cortexl.tri] = read_surf( fullfile(DirVal, 'lh.pial') );

%Draws the 3D figure of the brain with an opaque surface(ideal for Grid & Strips)
%Changing the "FaceAlpha" number will change the opacity
tripatch(cortexr,fig,[.5 .5 .5],'facealpha',.3);%Right Cortex
hold on
tripatch(cortexl,fig,[.5 .5 .5],'facealpha',.3);%Left Cortex

%Beautification
shading interp; lighting gouraud; material dull;

%Initializing Variables
ChanList=char(RAS_labels);%Channel List, a char array of all the channel names
ChannelName={};ChannelGrouping={};ChannelNumber=[];%Initializes empty cell arrays/chan_num array

%~~~~~~~~~Loops through every electrode contact in the channel list~~~~~~~~
for chan=1:length(ChanList)
    
    %currentChan = Char of current Elec Chan E.g. 'LOF01'
    currentChan=ChanList(chan,:);
    Num=[];
    
    %This FOR loop iterates through the
    for charsInChan=1:length(currentChan)
        %If the current iter char is a num
        if isstrprop(currentChan(charsInChan), 'digit')
            Num=[Num;charsInChan];%The indices with nums in the Label
        end
    end
    
    ChannelName{chan}=currentChan;
    
    %If the label is a micro, just skip it
    if strncmp('micro',currentChan,5)==1
        continue
    end
    
    %ChannelGrouping{chan}=currentChan(1:Num(1));
    ChannelGrouping{chan}=currentChan(1:Num(1)-1);
    ChannelNumber(chan)=str2num(currentChan(Num(1):end));
    %RASPerChan=RAS_coords(chan,:); ------UNUSED
    
end
%~~~~~~~~~Done Iterating through every Electrode Contact~~~~~~~~~~~~~~~~~~

%Creates Indices Matrix for each full Electrode, and applies the index
%notation to the full channel list so that each contact on the same
%electrode has the same index. (E.g. All LOF channels will have value of 3)
[IM,Label]=grp2idx(ChannelGrouping);
[MI,NonLRelecLabels]=grp2idx(ChanList(:,2:3));
COL=colormap(cool(max(MI)));
colormap(gray(255))
COL=colormap(cool(max(MI)));

% Change these to move the electrode contact itself
elecLR=0; %neg is left, pos is right
elecPA=0; %neg is posterior, pos is anterior
elecIS=0; %neg is Inferior, pos is Superior

%Change these to move the text for each contact
textLR=0; %neg is left, pos is right
textPA=0; %neg is posterior, pos is anterior
textIS=0; %neg is Inferior, pos is Superior

%This FOR loop iterates through each full Electrode and plots each contact
%in that electrode onto the 3D render of the brain
for electrode=1:length(Label)
    ElecNum=str2num(ChanList(IM==electrode,4:end)); %4:end MGH, 3:end BW
    scatter3(RAS_coords(IM==electrode,1)+elecLR,RAS_coords(IM==electrode,2)+elecPA,RAS_coords(IM==electrode,3)+elecIS,80,COL(MI(IM==electrode),:),'filled')
    hold on
    text(RAS_coords(IM==electrode,1)+textLR,RAS_coords(IM==electrode,2)+textPA,RAS_coords(IM==electrode,3)+textIS,ChanList(IM==electrode,:),'fontsize',10,'color',[0 0 0])
    
    % White Text--Acts on a certain subset of electrodes as determined by s1
    %     if sl>160
    %         ElecNum=str2num(CL(IM==sl,4:end)); %4:end MGH, 3:end BW
    %         scatter3(RAS_coords(IM==sl,1)+elecLR,RAS_coords(IM==sl,2)+elecPA,RAS_coords(IM==sl,3)+elecIS,80,COL(MI(IM==sl),:),'filled')
    %         hold on
    %         text(RAS_coords(IM==sl,1)+textLR,RAS_coords(IM==sl,2)+textPA,RAS_coords(IM==sl,3)+textIS,CL(IM==sl,:),'fontsize',10,'color',[1 1 1])
    %     end
end

axis vis3d%Makes the brain farther away in the figure view
daspect([1 1 1])
% title(['Patient: ', PatientName],'fontsize',24)
axis off

%~~~~~~~~~~~Pictures~~~~~~~~~~~~~~
V=[-90 0;90 0;180 90;-180 -90; 0 0];
for OL=1:length(V)
    view(V(OL,1),V(OL,2))
    lighting gouraud
    camlight('headlight')
    %
    %   shading interp; lighting gouraud; %material dull;
    %If you want to save this image:
    
    %  print(gcf,'-djpeg','-r600',['M:\Recons\MG103\',num2str(OL)])
    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% This is for the representations where the slices are along the axis of the electrode of interest.

%SKIP -- SKIP -- SKIP -- SKIP -- SKIP
%~~~~~~~~~A Bunch of Variables that are essentially the same thing, namely
%they are cell arrays with the names of the electrode contacts in it
ChanList=char(RAS_labels);

ChannelName={};ChannelGrouping={};ChannelNumber=[];
for chan=1:length(ChanList)
    currentChan=ChanList(chan,:);
    Num=[];
    for charsInChan=1:length(currentChan)
        if isempty(str2num(currentChan(charsInChan))) == 0 %IF the current iter char is a num
            Num=[Num;charsInChan]; % Constantly appended with the number of channels per electrode
        end
    end
    ChannelName{chan}=currentChan; %1 x N matrix of Channel Names
    ChannelGrouping{chan}=currentChan(1:Num(1)-1); %1 x N matrix of Channel Names
    ChannelNumber(chan)=str2num(currentChan(Num(1):end)); %1 x N Matrix of Channel Numbers per Electrode
    RASPerChan=RAS_coords(chan,:); %The RAS of the current iterations channel
end
[IM,Label]=grp2idx(ChannelGrouping); % IM is a list of all groups determined by grp2idx(useless?)
%~~~~~~~~

axis vis3d
daspect([1 1 1]) %No isometric viewing, all units on the axes are 1
title(['Patient: ',PatientName],'fontsize',24)
axis off

%~~~~~Takes pictures of the 3D Brain
V=[0 90;90 0;180 0];
for OL=1:1
    view(V(OL,1),V(OL,2))
    %     lighting none
    %     lighting gouraud
    camlight('headlight')
    %
    shading interp; lighting gouraud; material dull;
    %         If you want to save this image:
    %         print(gcf,'-djpeg','-r600',['C:\Users\COD Advanced Warfare\Documents\DARPA\CashLabECR\LocationViews\NoStim',Patient{IN},'View_',num2str(OL)])
end
%~~~~~~

%~~~~~~~Loading MRI File from nifti(.nii)
mri_file=fullfile(MRIDirectory,'mri.nii');
hdr_mr=load_nifti(mri_file);
vol = hdr_mr.vol;
%~~~~~~~~

%Defining a New RAS_coords1 for some reason...Change to RAS_coords
RAS_coords1=[RAS_coords(:,1) RAS_coords(:,2) RAS_coords(:,3)];

PLOT3D=1; % change this value if you don't want to do the glass brains
TargetChannels=[]; %this is if you want to label in red certain channel numbers

%~~~~~~~~Main Logic Loop from here to end of Cell~~~~~~~~~~~~~~
for ChanLabel=1:length(Label)%For each Depth/Grid/Strip
    
    ElecLabel=Label{ChanLabel};%This is the current working label
    
    %if  isempty(findstr('IOD',ElecLabel))==0 && isempty(findstr('REF',ElecLabel))==1 %&& isempty(findstr('FO',ElecLabel))==1
    clf
    figure1=gcf;
    set(gcf,'Position',[1 -29 1920 1004])
    Subplots=[1 4 7];
    for OL=1:3
        
        %~~~~~~~~~~~Plot Glass Brains~~~~~~~~~~~~~
        if PLOT3D==1
            subplot(3,3,Subplots(OL))
            patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[0.8 0.8 0.8],'facealpha',0.55,'edgecolor','none')
            hold on
            patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[0.8 0.8 0.8],'facealpha',0.55,'edgecolor','none')
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        %~~~~~~~~~Plotting the dots colored for the electrode depth~~~~
        [MI,NonLRelecLabels]=grp2idx(ChanList(:,2:3));%Indexification of the Electrodes by type
        COL=colormap(cool(max(MI)));%ColorMap of the last of the elec type???
        
        for electrode=1:length(Label)
            ElecNum=str2num(ChanList(IM==electrode,4:end));
            if isempty(strmatch(Label(electrode),ElecLabel))==0
                if PLOT3D==1
                    %                       scatter3(RAS_coords1(IM==sl,1),RAS_coords1(IM==sl,2),RAS_coords1(IM==sl,3),32,COL(MI(IM==sl),:),'filled')
                    scatter3(RAS_coords1(IM==electrode,1),RAS_coords1(IM==electrode,2),RAS_coords1(IM==electrode,3),32,'r','filled')
                    hold on
                end
                RASToPlot=RAS_coords1(IM==electrode,:);
                if PLOT3D==1
                    scatter3(RASToPlot(TargetChannels,1),RASToPlot(TargetChannels,2),RASToPlot(TargetChannels,3),32,'r','filled')
                end
            else
                if PLOT3D==1
                    scatter3(RAS_coords1(IM==electrode,1),RAS_coords1(IM==electrode,2),RAS_coords1(IM==electrode,3),16,[0.6 0.6 0.6],'filled')
                    hold on
                end
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        %~~~~~~~~~Adding the lines on 3D plot to see where you are~~~~~
        AngleAP=atand((RASToPlot(end,2)-RASToPlot(1,2))/(RASToPlot(end,1)+RASToPlot(1,1)));
        if PLOT3D==1
            hsp2= plot3(-100:100,RASToPlot(1,2)+zeros(length(-100:100),1),RASToPlot(1,3)+zeros(length(-100:100),1)-1,'-','linewidth',2,'color',[0.2 0.2 0.9]);
            rotate(hsp2,[0,0,1],AngleAP)
            hold on
            %                 plot3(-RASToPlot(1,1)+zeros(length(-100:100),1),RASToPlot(1,2)+zeros(length(-100:100),1),-100:100,'-','linewidth',2,'color',[0.9 0.2 0.2])
            plot3(RASToPlot(4,1)+zeros(length(-100:100),1),RASToPlot(4,2)+zeros(length(-100:100),1),-100:100,'-','linewidth',2,'color',[0.9 0.2 0.2])
            
            axis vis3d
            daspect([1 1 1])
            title([ElecLabel],'fontsize',14)
            axis off
            V=[0 90;90 0;0 0];
            
            view(V(OL,1),V(OL,2))
            lighting none
            lighting gouraud
            camlight('headlight')
            camlight('headlight')
            xlim([-80 80])
            ylim([-105 80])
            zlim([-80 100])
        end
        %~~~~~~~~~~~~~~~
        %end
        
        %%
        %Finding the MRI file
        mri_file=fullfile(MRIDirectory,'mri.nii');
        hdr_mr=load_nifti(mri_file);
        
        %transform RAS to vox
        vox2ras_mr=hdr_mr.vox2ras;
        ras2vox_mr=inv(vox2ras_mr);
        ras2vox=ras2vox_mr;
        ras2vox(1:3,4)=size(vol)'/2;
        Rdim=find(abs(ras2vox(:,1))==max(abs(ras2vox(:,1))));
        Adim=find(abs(ras2vox(:,2))==max(abs(ras2vox(:,2))));
        Sdim=find(abs(ras2vox(:,3))==max(abs(ras2vox(:,3))));
        
        pix=ras2vox*[RASToPlot, ones(size(RASToPlot,1),1)]';
        [numv,yloc]= grp2idx(round(pix(2,:)));
        [numv,dkd]= grp2idx(round(pix(3,:)));
        MidpointPerElectrodeX=nanmean(RASToPlot(:,1));
        MidpointPerElectrodeY=nanmean(RASToPlot(:,2));
        MidpointPerElectrodeZ=nanmean(RASToPlot(:,3));
        voxMid=ras2vox*[MidpointPerElectrodeX MidpointPerElectrodeY MidpointPerElectrodeZ 1]';
        voxMid=voxMid(1:3);
        
        AngleAP=atand((RASToPlot(end,2)-RASToPlot(1,2))/(RASToPlot(end,1)-RASToPlot(1,1)));
        AngleAPrad=atan((RASToPlot(end,2)-RASToPlot(1,2))/(RASToPlot(end,1)-RASToPlot(1,1)));
        Zmove=AngleAPrad*(128-str2num(char(dkd(1))));
        axes1 = axes('Parent',figure1,'Color',[0 0 0],'CLim',[0 255],...
            'Position',[0.2381 0.0320 0.6960 0.9388],...
            'DataAspectRatio',[1 1 1],'Clipping','on','Projection','orthographic'); %,'ClippingStyle','3dbox'
        xd=[];yd=[];zd=[];
        for IS=1:1
            %To create the slice to go through the MRI along the angle axis
            if length(dkd)>3
                hsp = surf(linspace(1,256,256),linspace(1,256,256),...
                    round( str2num(char(dkd(round(length(dkd)-length(dkd)/2)))))+zeros(256));
                %                 hsp = surf(linspace(1,256,256),linspace(1,256,256),...
                %                     round( str2num(char(dkd(1))))+zeros(256)+Zmove);
            else
                hsp = surf(linspace(1,256,256),linspace(1,256,256),...
                    round( str2num(char(dkd(1))))+zeros(256));
            end
            rotate(hsp,[1,0,0],-AngleAP) %Rotate along the angle axis in the x direction
            %               rotate(hsp,[0,0,1],-AngleAP) %Rotate along the angle axis in the z direction
            
            hold on
            xd(:,:,IS) = get(hsp,'XData');%hsp.XData;
            yd(:,:,IS) = get(hsp,'YData');%hsp.YData;
            zd(:,:,IS) = get(hsp,'ZData');%hsp.ZData;
            delete(hsp)
        end
        
        colormap(gray(255))
        subP=[2 3 5 6 8 9];
        hold on
        view(axes1,[90 90]);
        grid(axes1,'on');
        hold(axes1,'on');
        
        slice(vol,xd(:,:,1),yd(:,:,1),zd(:,:,1),axes1)
        hold on
        %     hsp3 = surf(linspace(1,256,256),linspace(1,256,256),...
        %         zeros(256));
        step=30;
        for si=1:length(dkd)
            chanNum=find(numv==si);
            PIn=pix(:,numv==si);
            for sim=1:length(chanNum)
                hold on
                %                 scatter3(PIn(2,sim),PIn(1,sim),PIn(3,sim)+10,32,[1 0 0], 'filled')
                scatter3(PIn(2,sim),PIn(1,sim),PIn(3,sim)+step,32,[1 0 0], 'filled')
            end
            for sim=1:length(chanNum)
                text(PIn(2,sim),PIn(1,sim),PIn(3,sim)+step+2,num2str(chanNum(sim)),'color','w','fontsize',9)
                hold on
            end
            shading flat
            axis image
            set(axes1,'xlim',[0 250],'ylim',[0 250])
            title([ElecLabel,' Channels: ',num2str(1),'-',num2str(length(numv))])
        end
        text(160,60,PIn(3,sim)+8,'R','color','w','fontsize',12,'backgroundcolor','k')
        text(160,180,PIn(3,sim)+8,'L','color','w','fontsize',12,'backgroundcolor','k')
        axis off
        %%
        %           pause
        colormap(gray(255))
        figure1.Color=[1 1 1];
        figure1.PaperPositionMode = 'auto';
        figure1.InvertHardcopy = 'off';
        print(gcf,'-dpng','-r400',[ImageDirectory,'ReconPerDepthSlice',ElecLabel])
        %            close
        
    end
end
%~~~~~~~~~End of Main Logic Loop~~~~~~~~~~~~~

% SKIP -- SKIP -- SKIP -- SKIP -- SKIP

%% This is for the single voxel representations where the slices are
% at the single voxel locations around the electrode of interest.
% CL=char(RAS_labels);

ChannelName={};ChannelGrouping={};ChannelNumber=[];
for chan=1:length(ChanList)
    
    currentChan=ChanList(chan,:);
    Num=[];
    
    for charsInChan=1:length(currentChan)
        %If the current iter char is a num
        if isstrprop(currentChan(charsInChan), 'digit')
            Num=[Num;charsInChan];%The indices with nums in the Label
        end
    end
    
    ChannelName{chan}=currentChan;
    ChannelGrouping{chan}=currentChan(1:Num(1)-1);
    ChannelNumber(chan)=str2num(currentChan(Num(1):end));
    RASPerChan(chan,:)=[RAS_coords(chan,1) RAS_coords(chan,2) RAS_coords(chan,3)];
end
[IM,Label]=grp2idx(ChannelGrouping);

TargetChannels=[];
PLOT3D=1; %Change if you don't want to run the glass brains

for CHa=1:length(Label) %this automatically runs through the labels for the depths
    % for CHa=2             %this automatically runs through the labels for the depths
    ElecLabel=Label{CHa};
    TargetElec=ElecLabel;
    if isempty(findstr('REF',ElecLabel))==1 %&& isempty(findstr('FO',Group))==1
        
        for ChannelNum=1:length(find(IM==CHa))
            clf
            
            figure1=gcf;
            set(gcf,'Position',[1 -3 1518 1000]);
            Subplots=[1 2 3];
            POS2=[[0.0500 0.7253 0.2134 0.2157];[0.7030 0.7253 0.2134 0.2157];[0.3530 0.7253 0.2134 0.2157];];
            AxesLim1=[[25 225];[25 225];[25 175]];
            AxesLim2=[[50 200];[25 200];[50 200]];
            for OL=1:3
                axes2 = axes('Parent',figure1,'Color',[0 0 0],'CLim',[100 255],...
                    'Position',[POS2(OL,:)],...
                    'DataAspectRatio',[1 1 1],'Clipping','on','Projection','orthographic');%,'ClippingStyle','3dbox');
                if PLOT3D==1
                    verticeslh1=verticeslh;
                    verticeslh1(:,1)=-verticeslh1(:,1);
                    patch('Faces',faceslh,'Vertices',verticeslh1,'FaceColor',[0.8 0.8 0.8],'facealpha',0.55,'edgecolor','none')
                    hold on
                    verticesrh1=verticesrh;
                    verticesrh1(:,1)=-verticesrh1(:,1);
                    patch('Faces',facesrh,'Vertices',verticesrh1,'FaceColor',[0.8 0.8 0.8],'facealpha',0.55,'edgecolor','none')
                    
                end
                
                COL=colormap(cool(max(MI)));
                for electrode=1:length(Label)
                    ElecNum=str2num(ChanList(IM==electrode,4:end));
                    if isempty(strmatch(Label(electrode),TargetElec))==0
                        if PLOT3D==1
                            scatter3(-RASPerChan(IM==electrode,1),RASPerChan(IM==electrode,2),RASPerChan(IM==electrode,3),20,COL(MI(IM==electrode),:),'filled')
                            hold on
                        end
                        RASToPlot=RASPerChan(IM==electrode,:);
                        RASToPlot=RASToPlot(ChannelNum,:);
                        if PLOT3D==1
                            scatter3(-RASToPlot(1,1),RASToPlot(1,2),RASToPlot(1,3),40,'r','filled')
                        end
                    else
                        if PLOT3D==1
                            scatter3(-RASPerChan(IM==electrode,1),RASPerChan(IM==electrode,2),RASPerChan(IM==electrode,3),16,[0.7 0.7 0.7],'filled')
                            hold on
                        end
                    end
                end
                
                if PLOT3D==1
                    hsp2= plot3(-100:100,RASToPlot(1,2)+zeros(length(-100:100),1),RASToPlot(1,3)+zeros(length(-100:100),1)-1,'-','linewidth',2,'color',[0.2 0.2 0.9]);
                    hold on
                    plot3(-RASToPlot(1,1)+zeros(length(-100:100),1),RASToPlot(1,2)+zeros(length(-100:100),1),-100:100,'-','linewidth',2,'color',[0.2 0.9 0.2])
                    axis vis3d
                    daspect([1 1 1])
                    title([TargetElec,' Channel: ',num2str(ChannelNum)],'fontsize',14)
                    axis off
                    V=[0 90;90 0;0 0];
                    view(V(OL,1),V(OL,2))
                    lighting none
                    lighting gouraud
                    camlight('headlight')
                    camlight('headlight')
                    xlim([-80 80])
                    ylim([-105 80])
                    zlim([-80 100])
                end
            end
            %
            
            %Finding the MRI file
            mri_file=fullfile(MRIDirectory,'mri.nii');
            hdr_mr=load_nifti(mri_file);
            vol = hdr_mr.vol;
            
            %transform RAS to vox
            vox2ras_mr=hdr_mr.vox2ras;
            ras2vox_mr=inv(vox2ras_mr);
            ras2vox=ras2vox_mr;
            ras2vox(1:3,4)=size(vol)'/2;
            Rdim=find(abs(ras2vox(:,1))==max(abs(ras2vox(:,1))));
            Adim=find(abs(ras2vox(:,2))==max(abs(ras2vox(:,2))));
            Sdim=find(abs(ras2vox(:,3))==max(abs(ras2vox(:,3))));
            
            RASToPlot2=[RASToPlot(:,1) RASToPlot(:,2) RASToPlot(:,3)];
            pix=ras2vox*[RASToPlot2, ones(size(RASToPlot,1),1)]';
            [numv,dkd]= grp2idx(round(pix(3,:)));
            MidpointPerElectrodeX=nanmean(RASToPlot(:,1));
            MidpointPerElectrodeY=nanmean(RASToPlot(:,2));
            MidpointPerElectrodeZ=nanmean(RASToPlot(:,3));
            voxMid=ras2vox*[MidpointPerElectrodeX MidpointPerElectrodeY MidpointPerElectrodeZ 1]';
            voxMid=voxMid(1:3);
            
            POS=[[0.001 0.0620 0.30 0.6];[0.303 0.0620 0.3 0.6];[0.603 0.0620 0.4 0.6]];
            
            %HERES HOW YOU CHANGE THE SIZES OF MRI STUFF
            AxesLim1=[[10 240];[10 250];[10 240]];%X-Direction
            AxesLim2=[[30 210];[30 210];[10 250]];%Y-Direction
            
            %this is to plot the different views from the MRI
            for VI=1:3
                axes1 = axes('Parent',figure1,'Color',[0 0 0],'CLim',[0 512],...
                    'Position',[POS(VI,:)],...
                    'DataAspectRatio',[1 1 1],'Clipping','on','Projection','orthographic');%,'ClippingStyle','3dbox');
                xd=[];yd=[];zd=[];
                
                for IS=1:3
                    if length(dkd)>3
                        hsp = surf(linspace(1,256,256),linspace(1,256,256),...
                            round( str2num(char(dkd(round(length(dkd)-length(dkd)/2)))))+zeros(256));
                    else
                        hsp = surf(linspace(1,256,256),linspace(1,256,256),...
                            round( str2num(char(dkd(1))))+zeros(256));
                    end
                    hold on
                    xd(:,:,IS) = get(hsp,'XData');%hsp.XData;
                    yd(:,:,IS) = get(hsp,'YData');%hsp.YData;
                    zd(:,:,IS) = get(hsp,'ZData');%hsp.ZData;
                    delete(hsp)
                end
                
                colormap(gray(255))
                subP=[2 3 5 6 8 9];
                hold on
                
                %     view(axes1,[0 0]);
                grid(axes1,'on');
                hold(axes1,'on');
                
                %This is the Recon Per Depth MRI picture Orientation
                %vol(:,:,:) corresponds to the volume in the Coronal,
                %Axial, and Sagittal Directions, respectively.
                %E.g.: vol(Coronal, Axial, Sagittal)
                if VI==1
                    imagesc(squeeze(vol(:,round(voxMid(2)),:)))
                    hold on
                    view(axes1,[90 90]);
                    text(20,45,'R','color','w','fontsize',12,'backgroundcolor','k')
                    text(20,190,'L','color','w','fontsize',12,'backgroundcolor','k')
                    
                elseif VI==2
                    imagesc(squeeze(vol(:,:,round(voxMid(3)))))
                    hold on
                    view(axes1,[90 90]);
                    %                     axis ij
                    text(20,40,'R','color','w','fontsize',12,'backgroundcolor','k')
                    text(20,190,'L','color','w','fontsize',12,'backgroundcolor','k')
                    
                elseif VI==3
                    imagesc(squeeze(vol(round(voxMid(1)),:,:)))
                    hold on
                    view(axes1,[-180 90]);
                    text(30,20,'A','color','w','fontsize',12,'backgroundcolor','k')
                    text(220,20,'P','color','w','fontsize',12,'backgroundcolor','k')
                    %                    figure1.Position=[1 -3 1518 820];
                end
                
                
                %To plot the right values relative to the channels and slices
                %This needs to be changed with the above if statements
                for si=1:length(dkd)
                    chanNum=find(numv==si);
                    PIn=pix(:,numv==si);
                    for sim=1:length(chanNum)
                        if VI==1
                            scatter(PIn(3,sim),PIn(1,sim),56,'r', 'filled');
                            hold on
                            text(PIn(3,sim),PIn(1,sim),num2str(ChannelNum),'color','w','fontsize',16)
                        elseif VI==2
                            title([TargetElec,' Channel: ',num2str(ChannelNum)],'fontsize',16)
                            scatter(PIn(2,sim),PIn(1,sim),56,'r', 'filled');
                            hold on
                            text(PIn(2,sim),PIn(1,sim),num2str(ChannelNum),'color','w','fontsize',16)
                        elseif VI==3
                            scatter(PIn(3,sim),PIn(2,sim),56,'r', 'filled');
                            hold on
                            text(PIn(3,sim),PIn(2,sim),num2str(ChannelNum),'color','w','fontsize',16)
                        end
                    end
                    
                    shading flat
                    set(axes1,'xlim',AxesLim1(VI,:),'ylim',AxesLim2(VI,:))
                    
                    
                    
                end
                axis off
                colormap(gray(255))
                caxis([0 300])
            end
            
            
            %This allows you to move the subplots in the figure(children)
            %around the layers, in this case, moving the rightmost MRI
            %behind the middle MRI
            k = get(gcf,'children');
            set(gcf,'children',k([2 1 3 4 5 6]))
            %
            %             pause;
            %             return;
            
            %~~~~~~~~~~~This section Seems to only work in Linux~~~~~~~~~~~
            figure1.Color=[1 1 1];
            figure1.PaperPositionMode = 'auto';
            figure1.InvertHardcopy = 'off';
            
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 32 20]);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            pause(2)
            print(gcf,'-dpng','-r400',fullfile( ImageDirectory, ['ReconPerDepth',TargetElec,'Channel',num2str(ChannelNum)]) )
            close
        end
    end
end


