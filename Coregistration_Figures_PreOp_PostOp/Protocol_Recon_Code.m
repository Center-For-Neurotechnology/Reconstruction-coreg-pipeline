%% Modular Reconstruction and Co-registration of Imaging from Implanted ECoG and SEEG Electrodes
% This version was completed on November 15, 2022
% 
% For questions, please post issues to the GitHub
% 

clear all
close

%Add the paths in this order only:
% Otherwise, functions in the fieldtrip toolbox which will make running this code
%difficult.
addpath(genpath('/MyPC/Documents/MATLAB/fieldtrip-20220202'))
addpath(genpath('/MyPC/bourbon-the-huckster/Documents/MATLAB/Recons/Coregistration_Figures_PreOp_PostOp'))

%Path to the recons directory. This should be filled with subject folders.
%Each subject folder will have the SurferOutput and RAS Excel 
wd = ['MyPC/Documents/Recons/'];

PatientName='sub-5o1r'; % name of the folder and id## in the working directory

% Setting this variable to 0 will create only the 2D visuals without
% required a SurferOutput from FreeSurfer
surferOutput = 1; % Set to 1 if you have run the Freesurfer Reconstruction

% This checks to see if you used BIDS formatting before continuing 
if isfile( fullfile( wd, PatientName, 'ses-postimp', 'ieeg', [PatientName '_ses-postimp_electrodes.tsv'] ) )
    BIDS = true;
elseif isfile( fullfile( wd, PatientName, [PatientName '_RAS.xlsx'] ) )
    BIDS = false;
else
    disp('We could not find the RAS excel spreadsheet in normal or BIDS format. Please check that all files are in the correct directories')
    return
end

if surferOutput == 0
    ImageDirectory = fullfile(wd , PatientName, 'Images'); %print files
    RASExcelDirectory = fullfile(wd, PatientName );
    MRIDirectory = fullfile(wd, PatientName ); %mri.nii...NOT nii.gz
    MRIFile=dir( fullfile(MRIDirectory, 'mri.nii'));
    [num,txt,all]=xlsread( fullfile( RASExcelDirectory, [PatientName '_RAS.xlsx'] ) );
    RAS_coords=[];
    RAS_coords=num;
    RAS_labels={};
    RAS_labels=txt(2:end,1);
elseif BIDS 
    DirVal = fullfile(wd , 'derivatives' , 'freesurfer', [PatientName,'_SurferOutput'], 'surf'); %the surf folder with the freesurfer output
    ImageDirectory = fullfile(wd , 'derivatives', 'ReconImages', PatientName, [PatientName,'_Images']);         %print files
    RASExcelDirectory = fullfile(wd, PatientName, 'ses-postimp','ieeg' );
    MRIDirectory = fullfile(wd , PatientName, 'ses-preimp', 'anat');
    MRIFile=dir( fullfile(MRIDirectory, '*.nii'));    
    RASInf = tsvread( fullfile(RASExcelDirectory, [PatientName '_ses-postimp_electrodes.tsv'] ) );
    RAS_coords=[];
    RAS_coords=[RASInf.x RASInf.y RASInf.z];
    RAS_labels={};
    RAS_labels=RASInf.name;
else
    DirVal = fullfile(wd , PatientName, [PatientName '_SurferOutput'], 'surf'); %the surf folder with the freesurfer output
    ImageDirectory = fullfile(wd , PatientName, 'Images');         %print files
    RASExcelDirectory = fullfile(wd, PatientName );
    MRIDirectory = fullfile(wd , PatientName, [PatientName '_SurferOutput'], 'mri'); %mri.nii...NOT nii.gz
    MRIFile=dir( fullfile(MRIDirectory, 'mri.nii'));
    [num,txt,all]=xlsread( fullfile( RASExcelDirectory, [PatientName '_RAS.xlsx'] ) );
    RAS_coords=[];
    RAS_coords=num;
    RAS_labels={};
    RAS_labels=txt(2:end,1);
end

% Check for Images Directory in the Working Directory, make it
% if it's not there
if ~isfolder(ImageDirectory)
    mkdir(ImageDirectory)
    disp('Made the Images folder in your Patient''s Working Directory!')    
end

display("Made it Through. Check the Workspace, that the files might match your expectations.");

% Getting the freesurfer pial and seeing where things are located, with the white matter in blue and pial surface.

% It's also handy for the blender-friendly stl file saving.
if isfile(fullfile(DirVal,'rh.pial.T1'))
    rightPialFile='rh.pial.T1';
else
    rightPialFile='rh.pial';
end

[verticesrh, facesrh] = freesurfer_read_surf( fullfile(DirVal,rightPialFile) );
% stlwrite( fullfile(DirVal,'right.stl'), facesrh, verticesrh)

if isfile(fullfile(DirVal,'lh.pial.T1'))
    leftPialFile='lh.pial.T1';
else 
    leftPialFile='lh.pial';
end
[verticeslh, faceslh] = freesurfer_read_surf( fullfile(DirVal,leftPialFile) );
% stlwrite( fullfile(DirVal,'left.stl'), faceslh, verticeslh)

[verticesrhwhite, facesrhwhite] = freesurfer_read_surf( fullfile(DirVal,'rh.white') );
% stlwrite( fullfile(DirVal,'rightWhite.stl'), facesrhwhite, verticesrhwhite)

[verticeslhwhite, faceslhwhite] = freesurfer_read_surf( fullfile(DirVal,'lh.white') );
% stlwrite([DirVal,'leftWhite.stl'], faceslhwhite, verticeslhwhite)

%Surfer Output Tester
patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[.7 .7 .7],'facealpha',.9,'edgecolor','none')
hold on
patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[.7 .7 .7],'facealpha',.9,'edgecolor','none')
camlight('headlight');

h=gcf;
[cortexr.vert, cortexr.tri] = read_surf( fullfile(DirVal, rightPialFile) );
[cortexl.vert, cortexl.tri] = read_surf( fullfile(DirVal, leftPialFile) );

tripatch(cortexr,h,[.5 .5 .5],'facealpha',1);
hold on
tripatch(cortexl,h,[.5 .5 .5],'facealpha',1);
camlight('headlight');
shading interp; lighting gouraud; material dull;


%% Plotting RAS electrodes on the Surfer Output
clf
close
% Create figure handle in fig
% figure('units','normalized','outerposition',[0 0 1 1])
% fig=gcf;
fig = figure('units','normalized','outerposition',[0 0 1 1]);

%Returns a matrix of vertex coordinates and face lists for the given .pial file
[cortexr.vert, cortexr.tri] = read_surf( fullfile(DirVal, rightPialFile) );
[cortexl.vert, cortexl.tri] = read_surf( fullfile(DirVal, leftPialFile) );

%Draws the 3D figure of the brain with an opaque surface(ideal for Grid & Strips)
%Changing the "FaceAlpha" number will change the opacity
tripatch(cortexr,fig,[.5 .5 .5],'facealpha',.2);%Right Cortex
hold on
tripatch(cortexl,fig,[.5 .5 .5],'facealpha',.2);%Left Cortex

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

    if isempty(Num)==1
        ChannelGrouping{chan}=currentChan;
        ChannelNumber(chan)=0;
    else
        ChannelGrouping{chan}=currentChan(1:Num(1)-1);
        ChannelNumber(chan)=str2num(currentChan(Num(1):end));
    end
    
    
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
    ElecNum=ChannelNumber(IM==electrode);
    if ~isnan(RAS_coords(IM==electrode,1))
        scatter3(RAS_coords(IM==electrode,1)+elecLR,RAS_coords(IM==electrode,2)+elecPA,RAS_coords(IM==electrode,3)+elecIS,80,COL(MI(IM==electrode),:),'filled')
        hold on
        text(RAS_coords(IM==electrode,1)+textLR,RAS_coords(IM==electrode,2)+textPA,RAS_coords(IM==electrode,3)+textIS,ChanList(IM==electrode,:),'fontsize',10,'color',[0 0 0],'Interpreter','none')
    end
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
    
    print(fig,'-dpng','-r400', fullfile(ImageDirectory,[PatientName '_ThreeDView_' num2str(OL)] ))
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
close

%% This is for the representations where the slices are along the axis of the electrode of interest.
clf
close
%~~~~~~~~~A Bunch of Variables that are essentially the same thing, namely
%they are cell arrays with the names of the electrode contacts in it
ChanList=char(RAS_labels);

ChannelName={};ChannelGrouping={};ChannelNumber=[];
for chan=1:length(ChanList)
    currentChan=ChanList(chan,:);
    Num=[];
    for charsInChan=1:length(currentChan)
        if isstrprop(currentChan(charsInChan), 'digit') %IF the current iter char is a num
            Num=[Num;charsInChan]; % Constantly appended with the number of channels per electrode
        end
    end
    ChannelName{chan}=currentChan; %1 x N matrix of Channel Names
    %If the label is a micro, just skip it
    %     if strncmp('micro',currentChan,5)==1
    %         continue
    %     end
    if isempty(Num)==1
        ChannelGrouping{chan}=currentChan;
        ChannelNumber(chan)=0;
    else
        ChannelGrouping{chan}=currentChan(1:Num(1)-1);
        ChannelNumber(chan)=str2num(currentChan(Num(1):end));
    end
    RASPerChan=RAS_coords(chan,:); %The RAS of the current iterations channel
end
[IM,Label]=grp2idx(ChannelGrouping); % IM is a list of all groups determined by grp2idx(useless?)
%~~~~~~~~
fig = figure('units','normalized','outerposition',[0 0 1 1]);

axis vis3d
daspect([1 1 1]) %No isometric viewing, all units on the axes are 1
axis off

%~~~~~~

% This loads the MRI and then rotates it and can interpolate it to let the
% reslicing happen as it is easier to reslice along the length of the
% electrode if the input MRI has equal slices. The color scale is also adjusted.

% HOWEVER, use interpolation with caution as it could shift the center of the volume.
% If you want to adjust for coordinate systems (ACPC) use ft_volumerealign which
% is detailed here:
% https://flashsherlock.gitee.io/2021/01/12/ieeg-localization-of-intracranial-electrodes/
% and here:
% https://github.com/fieldtrip/fieldtrip/blob/master/ft_volumerealign.m

mri_file=fullfile(MRIDirectory,MRIFile(1).name);
mri = ft_read_mri(mri_file);
cfg  = [];
if length(unique(mri.dim))==1 %this checks if all the dimensions are the same
    cfg.method='flip'; %Flips orientation to be the same but no interpolation
else
    cfg.method='linear';  %Interpolates
end
cfg.downsample=1;
mrirs2 = ft_volumereslice(cfg,mri);
Bvol2=mrirs2.anatomy;
vol=255*(Bvol2/max(max(max(Bvol2))));
vox2ras_mr=mrirs2.transform;
ras2vox_mr=inv(vox2ras_mr);
% Slices of MRI relative to the whole electrode depth

%Defining a New RAS_coords1 for some reason...Change to RAS_coords
RAS_coords1=[RAS_coords(:,1) RAS_coords(:,2) RAS_coords(:,3)];

PLOT3D=surferOutput; % change this value if you don't want to do the glass brains
TargetChannels=[]; %this is if you want to label in red certain channel numbers

%~~~~~~~~Main Logic Loop from here to end of Cell~~~~~~~~~~~~~~
for ChanLabel=1:length(Label)%For each Depth/Grid/Strip
    
    ElecLabel=Label{ChanLabel};%This is the current working label
    
    Subplots=[1 4 7];
    for OL=1:1
        
        subplot(1,2,1)
        %~~~~~~~~~~~Plot Glass Brains~~~~~~~~~~~~~
        if PLOT3D==1
            %             subplot(3,3,Subplots(OL))
            patch('Faces',faceslh,'Vertices',verticeslh,'FaceColor',[0.8 0.8 0.8],'facealpha',0.2,'edgecolor','none')
            hold on
            patch('Faces',facesrh,'Vertices',verticesrh,'FaceColor',[0.8 0.8 0.8],'facealpha',0.2,'edgecolor','none')
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
            plot3(RASToPlot(ceil(size(RASToPlot,1)/2),1)+zeros(length(-100:100),1),RASToPlot(ceil(size(RASToPlot,1)/2),2)+zeros(length(-100:100),1),-100:100,'-','linewidth',2,'color',[0.9 0.2 0.2])
            
            axis vis3d
            daspect([1 1 1])
            title([ElecLabel],'fontsize',14,'Interpreter','none')
            axis off
            V=[0 90;90 0;0 0];
            
            view(-180,0)
            lighting none
            lighting gouraud
            camlight('headlight')
            camlight('headlight')
            xlim([-100 100])
            ylim([-100 100])
            zlim([-100 100])
        end
        axis off
        
        
        ras2vox=ras2vox_mr;
        ras2vox(1:3,4)=size(vol)'/2;
        
        pix=ras2vox*[RASToPlot, ones(size(RASToPlot,1),1)]';
        [numv,dkd]= grp2idx(round(pix(3,:)));
        
        axes1 = axes('Parent',fig,'Color',[0 0 0],'CLim',[0 255],...
            'Position',[0.4 -.02 0.6960 0.9388],...
            'DataAspectRatio',[1 1 1],'Clipping','on','Projection','orthographic'); %,'ClippingStyle','3dbox'
        
        colormap(gray(255))
        hold on
        RAS1= RASToPlot(1,:);
        RAS2= RASToPlot(end,:);
        
        ras2vox_mr=inv(vox2ras_mr);
        
        ras2vox=ras2vox_mr;
        ras2vox(1:3,4)=size(vol)'/2;
        vox1=ras2vox*[RAS1 1]';
        vox1=vox1(1:3);
        vox2=ras2vox*[RAS2 1]';
        vox2=vox2(1:3);
        
        %find normal of plane in vox
        Vslice=vox1-vox2;
        Vsup=ras2vox*[0 0 1 0]';
        Vsup=Vsup(1:3);
        Nplane=cross(Vslice,Vsup);
        Nplane=Nplane/norm(Nplane);
        
        %find rotation vector and angle
        Nz=[0 0 1];
        Nrotate=cross(Nplane,Nz);
        Nrotate=Nrotate/norm(Nplane);
        theta=acosd(dot(Nplane,Nz))/(norm(Nrotate)*norm(Nz));
        
        hsp=surf(ones(size(vol,1)*2,size(vol,2)*2)*vox2(3));
        
        rotate(hsp,Nrotate,theta,vox2);
        xd=get(hsp,'XData');
        yd=get(hsp,'YData');
        zd=get(hsp,'ZData');
        delete(hsp);
        
        h=slice(vol,yd,xd,zd);
        set(h,'edgecolor','none');
        new_slice=get(h,'cdata');
        [xarea,yarea]=find(~isnan(new_slice));
        new_slice=new_slice(min(xarea):max(xarea),min(yarea):max(yarea));
        
        hold on
        step=-5;
        for si=1:size(RASToPlot,1)
            chanNum=find(numv==si);
            PIn=pix(:,numv==si);
            for sim=1:length(chanNum)
                hold on
                scatter3(PIn(2,sim)-2,PIn(1,sim)+0,PIn(3,sim),100,[1 0 0], 'filled')
            end
            for sim=1:length(chanNum)
                text(PIn(2,sim)+step,PIn(1,sim),PIn(3,sim),num2str(chanNum(sim)),'color','w','fontsize',16)
                hold on
            end
            shading flat
            axis image
            set(axes1,'xlim',[0 250],'ylim',[0 250])
            
            view(-90,0)
            %             xlim([-80 80])
            ylim([36 220])
            zlim([36 220])
            
            if ~isnan(nanmean(new_slice(:)))
                caxis([ 0 nanmean(new_slice(:))+3*nanstd(new_slice(:))])
            else
                caxis([ 0 220])
            end
            title([ElecLabel,' Channels: ',num2str(1),'-',num2str(length(numv))],'fontsize',20,'Interpreter','none')
        end
        text(10,60,128,'L','color','w','fontsize',24,'backgroundcolor','k')
        text(10,210,128,'R','color','w','fontsize',24,'backgroundcolor','k')
        %         axis off
        set(gca,'xticklabel',{[]},'yticklabel',{[]},'zticklabel',{[]})
        
        colormap(gray(255))

        print(gcf,'-dpng','-r400',fullfile( ImageDirectory, ['ReconPerDepthSlice',ElecLabel]))
        clf
    end
end

close all

%% This is for the single voxel representations where the electrodes are
% at the single voxel locations  with the MRI slices around the electrode of interest.

% This loads the MRI and then rotates it without interpolating to use the
% raw data. The color scale is also adjusted.

mri_file=fullfile(MRIDirectory,MRIFile(1).name);
mri = ft_read_mri(mri_file);
cfg  = [];
if length(unique(mri.dim))==1 %this checks if all the dimensions are the same
    cfg.method='flip'; %Flips orientation to be the same but no interpolation
else
    cfg.method='linear';  %Interpolates
end
cfg.downsample=1;
mrirs = ft_volumereslice(cfg,mri);
Bvol=mrirs.anatomy;
vol2=256*Bvol/max(max(max(Bvol)));
vol2(isnan(vol2))=0;

vox2ras_mr2=mrirs.transform;
ras2vox_mr2=inv(vox2ras_mr2);

ChanList=char(RAS_labels);%Channel List, a char array of all the channel names
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
    if isempty(Num)==1
        ChannelGrouping{chan}=currentChan;
        ChannelNumber(chan)=0;
    else
        ChannelGrouping{chan}=currentChan(1:Num(1)-1);
        ChannelNumber(chan)=str2num(currentChan(Num(1):end));
    end
    RASPerChan(chan,:)=[RAS_coords(chan,1) RAS_coords(chan,2) RAS_coords(chan,3)];
end
[IM,Label]=grp2idx(ChannelGrouping);

TargetChannels=[];
PLOT3D=surferOutput; %Change if you don't want to run the glass brains

for CHa=1:length(Label) %this automatically runs through the labels for the depths
    % for CHa=2             %this automatically runs through the labels for the depths
    ElecLabel=Label{CHa};
    TargetElec=ElecLabel;
    if isempty(findstr('REF',ElecLabel))==1 %&& isempty(findstr('FO',Group))==1
        
        for ChannelNum=1:length(find(IM==CHa))
            %         for ChannelNum=1:1
            
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            
%             set(gcf,'Position',[1 -3 1518 1000]);
            Subplots=[1 2 3];
            POS2=[[0.0500 0.7253 0.2134 0.2157];[0.7030 0.7253 0.2134 0.2157];[0.3530 0.7253 0.2134 0.2157];];
            AxesLim1=[[25 225];[25 225];[25 175]];
            AxesLim2=[[50 200];[25 200];[50 200]];
            for OL=1:3
               
                if PLOT3D==1
                     axes2 = axes('Parent',fig,'Color',[0 0 0],'CLim',[100 255],...
                    'Position',[POS2(OL,:)],...
                    'DataAspectRatio',[1 1 1],'Clipping','on','Projection','orthographic');%,'ClippingStyle','3dbox');
                    verticeslh1=verticeslh;
                    verticeslh1(:,1)=-verticeslh1(:,1);
                    patch('Faces',faceslh,'Vertices',verticeslh1,'FaceColor',[0.8 0.8 0.8],'facealpha',0.55,'edgecolor','none')
                    hold on
                    verticesrh1=verticesrh;
                    verticesrh1(:,1)=-verticesrh1(:,1);
                    patch('Faces',facesrh,'Vertices',verticesrh1,'FaceColor',[0.8 0.8 0.8],'facealpha',0.55,'edgecolor','none')
                    
                end
                
                [MI,NonLRelecLabels]=grp2idx(ChanList(:,2:3));
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
                    V=[0 0;270 0;0 90];
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
            
            %transform RAS to vox
            vox2ras_mr=mrirs.transform;
            ras2vox_mr=inv(vox2ras_mr);
            ras2vox=ras2vox_mr;
            ras2vox(1:3,4)=size(vol2)'/2;
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
            
            POS=[[0.003 0.0620 0.3 0.6];[0.303 0.0620 0.3 0.6];[0.603 0.0620 0.4 0.6]];
            
            %HERE'S HOW YOU CHANGE THE SIZES OF MRI STUFF
            Coronal=squeeze(vol2(:,round(voxMid(2)),:));
            Sagittal=squeeze(vol2(round(voxMid(1)),:,:));
            Horizontal=squeeze(vol2(:,:,round(voxMid(3))));
            
            CoronalRange=squeeze(nanmean(vol2,2));
            SagittalRange=squeeze(nanmean(vol2,1));
            HorizontalRange=squeeze(nanmean(vol2,3));
            
            for smw=1:3
                if smw==1
                    rngAx=CoronalRange;
                    ImageRangeThresh=nanmedian(CoronalRange(:))-nanmedian(CoronalRange(:))*.6;
                elseif smw==2
                    rngAx=HorizontalRange;
                    ImageRangeThresh=nanmedian(HorizontalRange(:))-nanmedian(HorizontalRange(:))*.6;
                elseif smw==3
                    rngAx=SagittalRange;
                    ImageRangeThresh=nanmedian(SagittalRange(:))-nanmedian(SagittalRange(:))*.6;
                end
                
                IndZero=find(nanmean(rngAx)<ImageRangeThresh);
                first = IndZero(diff(IndZero)>1);
                second = IndZero(find(diff(IndZero)>1)+1);
                YRng=[first(1) second(1)];
                
                IndZero=find(nanmean(rngAx')<ImageRangeThresh);
                first = IndZero(diff(IndZero)>1);
                second = IndZero(find(diff(IndZero)>1)+1);
                XRng=[first(1) second(1)];
                
                if isempty(YRng)==1
                    YRng=[1 size(rngAx,2)];
                end
                if isempty(XRng)==1
                    XRng=[1 size(rngAx,1)];
                end
                
                AxesLim1(smw,:)=YRng;
                AxesLim2(smw,:)=XRng;

            end
            
            %this is to plot the different views from the MRI
            for VI=1:3
                axes1 = axes('Parent',fig,'Color',[0 0 0],'CLim',[0 512],...
                    'Position',[POS(VI,:)],...
                    'DataAspectRatio',[1 1 1],'Clipping','on','Projection','orthographic');%,'ClippingStyle','3dbox');
                colormap(gray(255))
                subP=[2 3 5 6 8 9];
                hold on
                
                %     view(axes1,[0 0]);
                grid(axes1,'on');
                hold(axes1,'on');
                
                %This is the Recon Per Depth MRI picture Orientation
                %vol2(:,:,:) corresponds to the volume in the Coronal,
                %Axial, and Sagittal Directions, respectively. 
                %E.g.: vol(Coronal, Axial, Sagittal)
                if VI==1
                    imagesc(Coronal)
                    hold on
                    view(axes1,[-90 90]);
                    text(AxesLim1(VI,2)*0.95,AxesLim2(VI,1)+10,'L','color','w','fontsize',12,'backgroundcolor','k')
                    text(AxesLim1(VI,2)*0.95,AxesLim2(VI,2)-10,'R','color','w','fontsize',12,'backgroundcolor','k')
                    shading flat
                    set(axes1,'xlim',AxesLim1(VI,:),'ylim',AxesLim2(VI,:))
                                        daspect([fliplr(size(Coronal)/max(size(Coronal))) 1])
                elseif VI==2
                    imagesc(Horizontal)
                    hold on
                    view(axes1,[270 90]);
                    %                  axis ij
                    text(AxesLim1(VI,2)*0.95,AxesLim2(VI,1)+10,'L','color','w','fontsize',12,'backgroundcolor','k')
                    text(AxesLim1(VI,2)*0.95,AxesLim2(VI,2)-10,'R','color','w','fontsize',12,'backgroundcolor','k')
                    shading flat
                    set(axes1,'xlim',AxesLim1(VI,:),'ylim',AxesLim2(VI,:))
                    %                     daspect([fliplr(size(Horizontal)/max(size(Horizontal))) 1])
                elseif VI==3
                    imagesc(Sagittal)
                    hold on
                    
                    view(axes1,[-90 90]);
                    text(AxesLim1(VI,2)*0.95,AxesLim2(VI,1)+10,'P','color','w','fontsize',12,'backgroundcolor','k')
                    text(AxesLim1(VI,2)*0.95,AxesLim2(VI,2)-10,'A','color','w','fontsize',12,'backgroundcolor','k')
                    shading flat
                    set(axes1,'xlim',AxesLim1(VI,:),'ylim',AxesLim2(VI,:))
%                     set(axes1, 'xlim', [1 300], 'ylim', AxesLim2(VI,:))
                    %                                         daspect([fliplr(size(Sagittal)/max(size(Sagittal))) 1])
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
                end
                
                axis off
                colormap(gray(255))
                VolCheck=vol2(:,:,round(size(vol2,3)/2));
                brightness = mean(VolCheck(:))+3*std(VolCheck(:));
                caxis([0 brightness])
            end
            
            pause(2)
  
            print(fig,'-dpng','-r400',fullfile( ImageDirectory, ['ReconPerDepth',TargetElec,'Channel',num2str(ChannelNum)]) )
            close
        end
    end
end
