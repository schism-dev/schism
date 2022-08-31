function []=SCHISM_SLAB2(base,varname,ivs,s_or_z,lev_or_zcor,ispher_nowrap,stacks,nspool,test)
% Authors: Joseph Zhang
% Date: Nov 2021
% Matlab function to visualize horizontal slabs at either a fix z level
% or at a fix sigma level from scribed outputs. Works for node-centered variables only!
% For a fixed z level, nan is used for above surface/below bottom.
% Works for mixed grid.
% Always requires out2d*.nc,  zCoordinates*.nc and get_global_info2.m
% SCHISM_SLAB2(base,varname,ivs,s_or_z,lev_or_zcor,ispher_nowrap,stacks,nspool,test)
% Inputs: 
%         base: base directory where base/outputs/ contains nc outputs *.nc
%         varname = var names like 'salinity' etc (node based only)
%         ivs: 1 for scalar, 2 for vector
%         s_or_z = 'S' (along a fixed sigma level) or 'Z' (along a fix z level).
%                  This is not used for 2D variables.
%         lev_or_zcor = level index (1 to nvrt) if s_or_z='S'; z-coordinate value 
%                       (z<0 is below MSL) if s_or_z='Z'. 
%                       This is not used for 2D variables.
%         ispher_nowrap: for spherical grid only. 1: remove wrap-around elem's across dateline in scalar display. 0: original
%         stacks: array of stack numbers (e.g. [2 4 5]) in the output file names (related to time)
%         nspool: sub-sampling frequency within each stack (e.g. '1' - include all)
%         test: 'y' (only plot out 1st frame for test); 'n' (plot all frames)
%         May need to adjust some parameters inside (e.g. caxis) to get right appearance of images
% Outputs: images and slab.avi

close all; 
%scrsz(1:2)=1; scrsz(3)=width (pixels), scrsz(4)=height
scrsz = get(0,'ScreenSize'); 
%4 parameters of position: left bottom_coord width height
figure('Position',[1 scrsz(4)*0.2 scrsz(3)/2 scrsz(4)*0.7]); 

% Read the variable and the vertical grid
delete('slab.avi');
vidObj = VideoWriter('slab.avi');
vidObj.FrameRate = 30;  % Default 30; smaller ->slower
%vidObj.Quality = 50;    % Default 75
open(vidObj);

if(strcmp(test,'y'))
  stacks2=stacks(1); 
else
  stacks2=stacks; 
end

%Get basic info
[file_nm,i23d,ne,np,nvrt,nm,xy00,vdry_e,h0,dp,kbp00,vid_eta]=get_global_info2(base,varname,stacks);

%xy00(np,:)
%nm(:,end)
%vzcor

%Mark all nodes in wrap-around elements for scalar display later
if(ispher_nowrap==1)
  %If a node is on dateline, make it 180 deg
  xtmp=xy00(:,1);
  xtmp(find(xtmp==-180))=180;
  xy00(:,1)=xtmp;
  icolor_nd(1:np)=0;
  for i=1:ne
    lon_min=1.e10; lon_max=-1.e10;
    if(isnan(nm(4,i))) 
      i34=3;
    else
      i34=4;
    end
    lon_min=min(xy00(nm(1:i34,i),1));
    lon_max=max(xy00(nm(1:i34,i),1));
    if(lon_min<-100 && lon_max>100)
      icolor_nd(nm(1:i34,i))=1;
    end
  end %for i
end %ispher_nowrap==1

%For plotting vector scales
xmax=max(xy00(:,1)); 
xmin=min(xy00(:,1)); 
ymax=max(xy00(:,2));
ymin=min(xy00(:,2));

%Prep output array
if(strcmp(s_or_z,'S') || i23d==1)
  out5=zeros(ivs,np);
else
  out5=zeros(ivs,nvrt,np);
  %zcor=zeros(nvrt,np);
end

istep=0;
for day=stacks2
  istep=istep+1;

  fname=[base,'/outputs/' file_nm '_' num2str(day) '.nc'];

  disp(['doing ' fname]);
  ncid1 = netcdf.open(fname,'NC_NOWRITE');
  vid1=netcdf.inqVarID(ncid1,varname);
  
  if(ivs==2)
    fname=[base,'/outputs/' file_nm(1:end-1) 'Y_' num2str(day) '.nc'];
    ncid2 = netcdf.open(fname,'NC_NOWRITE');
    vid2=netcdf.inqVarID(ncid2,[varname(1:end-1) 'Y']);
  end

  vid=netcdf.inqVarID(ncid1,'time');
  timeout0=double(netcdf.getVar(ncid1, vid)); %sec
  nrec=length(timeout0);

  if(strcmp(test,'y'))
    it2=1;
  else
    it2=nrec;
  end

  %Open up out2d and zCoordinates
  fname_2d=[base,'/outputs/out2d_' num2str(day) '.nc'];
  ncid_2d = netcdf.open(fname_2d,'NC_NOWRITE');

  fname_z=[base,'/outputs/zCoordinates_' num2str(day) '.nc'];
  ncid_z = netcdf.open(fname_z,'NC_NOWRITE');
  vid_z=netcdf.inqVarID(ncid_z,'zCoordinates');

  for it=1:nspool:it2;
    timeout=timeout0(it);
    time_d=fix(timeout/86400); %days
    time_h=fix((timeout-time_d*86400)/3600);
    time_m=fix((timeout-time_d*86400-time_h*3600)/60);
    time_s=timeout-time_d*86400-time_h*3600-time_m*60;

    %out5(ivs,np) if 2D or 'S'; otherwise out5(ivs,nvrt,np)
      eta2=double(netcdf.getVar(ncid_2d,vid_eta,[0 it-1],[np 1]));
      if(i23d==1) %2D
        out5(1,:)=double(netcdf.getVar(ncid1,vid1,[0 it-1],[np 1])); %(np)
        if(ivs==2)
          out5(2,:)=double(netcdf.getVar(ncid2,vid2,[0 it-1],[np 1]));
        end
      elseif(strcmp(s_or_z,'S'))
        out5(1,:)=double(netcdf.getVar(ncid1,vid1,[lev_or_zcor-1 0 it-1],[1 np 1])); %(np)
        if(ivs==2)
          out5(2,:)=double(netcdf.getVar(ncid2,vid2,[lev_or_zcor-1 0 it-1],[1 np 1]));
        end
      elseif(strcmp(s_or_z,'Z'))
        out5(1,:,:)=double(netcdf.getVar(ncid1,vid1,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
        if(ivs==2)
          out5(2,:,:)=double(netcdf.getVar(ncid2,vid2,[0 0 it-1],[nvrt np 1]));
        end
        zcor=double(netcdf.getVar(ncid_z,vid_z,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
      else
        error('Unknown s_or_z');
      end %i23d

%      elseif(ivs==2)
%        if(i23d==1) %2D
%          out5=double(netcdf.getVar(ncid1,vid1,[0 0 it-1],[2 np 1])); %(ivs,np)
%        elseif(strcmp(s_or_z,'S'))
%          out5=double(netcdf.getVar(ncid1,vid1,[0 lev_or_zcor-1 0 it-1],[2 1 np 1])); %(ivs,np)
%        elseif(strcmp(s_or_z,'Z'))
%          out5=double(netcdf.getVar(ncid1,vid1,[0 0 0 it-1],[2 nvrt np 1])); %(ivs,nvrt,np)
%          zcor=double(netcdf.getVar(ncid1,vzcor,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
%        else
%          error('Unknown s_or_z');
%        end
%      else
%        error('Unknown ivs');
%      end %ivs

    %Construct output uout(ivs,1:np) 
    if(i23d==1 || strcmp(s_or_z,'S'))
      uout=out5;
      %Mask dry nodes
      indx=find(eta2+dp<=h0);
      uout(:,indx)=nan;
    else
      %Interp in vertical @ wet spots only
      uout=NaN(ivs,np);
      for ii=1:np
        if(eta2(ii)+dp(ii)>h0)
%          indx=find(~isnan(zcor(:,ii)) & abs(zcor(:,ii))<1.e15);
%          indx2=indx(find(indx>=kbp00(ii)));
          for j=1:ivs
            uout(j,ii)=interp1(zcor(kbp00(ii):nvrt,ii),out5(j,kbp00(ii):nvrt,ii),lev_or_zcor,'linear',nan);
          end %j
        end %if
      end %ii
    end

    %Define junk
    uout(find(abs(uout)>1.e20))=nan;

    % plot uout
    % set axes
    v2=[xmin xmax ymin ymax]; %axis;
    % Write time stamp info
    loc_info_x=(v2(2)+v2(1))/2;
    loc_info_y=v2(4)*0.97+v2(3)*0.03;

    if(ivs==1) %scalar
      uout_p=uout(1,:);
      %Remove elem across dateline
      if(ispher_nowrap==1)
        uout_p(find(icolor_nd==1))=nan;
      end %ispher_nowrap

      h1=patch('Faces',nm','Vertices',xy00,'FaceVertexCData',uout_p','FaceColor','interp','EdgeColor','none');
      colormap(jet(40));
      % Set colormap range
      caxis([0 30]); colorbar;
    else %vector
      loc_scale_x=v2(2)*0.3+v2(1)*0.7;
      loc_scale_y=-v2(4)*0.02+v2(3)*1.02;
      x_aug=[xy00(:,1)' loc_scale_x];
      y_aug=[xy00(:,2)' loc_scale_y];
      quiver(x_aug,y_aug,[uout(1,:) 1],[uout(2,:) 0]);
      text(loc_scale_x,loc_scale_y,'1 m/s');

      %Alternatively, plot vector magnitude
%      h1=patch('Faces',nm','Vertices',xy00,'FaceVertexCData',sqrt(uout(1,:).^2+uout(2,:).^2)','FaceColor','interp','EdgeColor','none');
%      colormap(jet(40));
%      caxis([0 1]); colorbar;

    end %ivs

    axis([xmin xmax ymin ymax]);
    %axis([3.2e5 3.7e5 2.5e5 3.1e5]);
    text(loc_info_x,loc_info_y,{'Time (DD:HH:MM:SS)'; num2str([time_d time_h time_m time_s])});

    %axis off;

    % Add image to avi file
    set(gcf,'Color',[1 1 1]);
    set(gcf,'nextplot','replacechildren');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    if(day ~= stacks2(end) || it ~=it2)
      clf; %clear figure to avoid overlay
    end

  end %it
  netcdf.close(ncid1);
  netcdf.close(ncid_2d);
  netcdf.close(ncid_z);
  if(ivs==2); netcdf.close(ncid2); end;
end %for day=stacks2
close(vidObj);
