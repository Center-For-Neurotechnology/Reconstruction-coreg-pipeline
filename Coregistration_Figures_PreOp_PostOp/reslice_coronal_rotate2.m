%---------------------------------------------------
% FUNCTION: function slice = reslice_coronal(fname,RAS1,RAS2)
% INPUTS:   
%           fname = full path & filename of MR volume (nifti) of interest
%           
%           RAS1 = coordinate of most lateral (medial) electrode of single
%                  depth array
%           RAS2 = coordiante of most medial (lateral) electrode of single
%                  depth array
%                   
% OUTPUT:   Single semi-coronal slice aligned to array axis
%          
%---------------------------------------------------
% written by: Andrew Dykstra, MIT, 2009
% last edited: Alex Chan, 09/23/2009
% messed with by Angelique Paulk 9/14/2022
%--------------------------------------------------- 

% function [new_slice ras2vox] = reslice_coronal_rotate2(mri_file,ct_file,RAS1,RAS2)

ras2vox_mr=inv(vox2ras_mr);

ras2vox=ras2vox_mr;
        ras2vox(1:3,4)=size(vol)'/2;
% ras2vox(1:3,4)=ct2mr_vox*abs(ras2vox_ct(1:3,4)).*(size(hdr_mr.vol)'./abs(ct2mr_vox*size(hdr_ct.vol)'));
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

% hsp=surf(1:size(vol,1)*2,(1:size(vol,2)*2),ones(size(vol,1)*2,size(vol,2)*2)*vox2(3));
hsp=surf(ones(size(vol,1)*2,size(vol,2)*2)*vox2(3));

rotate(hsp,Nrotate,theta,vox2);
xd=get(hsp,'XData');
yd=get(hsp,'YData');
zd=get(hsp,'ZData');
delete(hsp);

h=slice(vol,yd,xd,zd);
set(h,'edgecolor','none');
%hold on;
%h=slice(vol,[20 size(vol,1)-20],[20],[20]);
%set(h,'edgecolor','none');
new_slice=get(h,'cdata');
[xarea,yarea]=find(~isnan(new_slice));
new_slice=new_slice(min(xarea):max(xarea),min(yarea):max(yarea));
