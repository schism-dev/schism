function []=SCHISM_SLAB(icomb,base,varname,s_or_z,lev_or_zcor,stacks,nspool,test)
% Authors: Joseph Zhang
% Date: Feb 2019
% Matlab function to visualize horizontal slabs at either a fix z level
% or at a fix sigma level. Works for node-centered variables only!
% For a fixed z level, nan is used for above surface/below bottom.
% Works for mixed grid.
% Requires get_global_info.m (in this dir)
% SCHISM_SLAB(icomb,base,varname,s_or_z,lev_or_zcor,stacks,nspool,test)
% Inputs: 
%         icomb: work with uncombined (0) or combined (1) nc 
%         base: base directory where base/outputs/ contains nc outputs schout_*.nc
%         varname = var names like 'salt' etc (node based only)
%         s_or_z = 'S' (along a fixed sigma level) or 'Z' (along a fix z level).
%                  This is not used for 2D variables.
%         lev_or_zcor = level index (1 to nvrt) if s_or_z='S'; z-coordinate value 
%                       (z<0 is below MSL) if s_or_z='Z'. 
%                       This is not used for 2D variables.
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
[ne,np,nvrt,nm,xy00,i23d,ivs,vzcor,vid5,vdry_e,h0,dp,kbp00,vid_eta,nproc,np_lcl,ne_lcl,ns_lcl,iplg,ielg,iegl_rank]=get_global_info(icomb,base,varname,stacks);

%xy00(np,:)
%nm(:,end)
%vzcor

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

  if(icomb==0) 
    fname=[base,'/outputs/' 'schout_0000_' num2str(day) '.nc'];
  else
    fname=[base,'/outputs/' 'schout_' num2str(day) '.nc'];
  end
  disp(['doing ' fname]);
  ncid0 = netcdf.open(fname,'NC_NOWRITE');

  vid=netcdf.inqVarID(ncid0,'time');
  timeout0=double(netcdf.getVar(ncid0, vid)); %sec
  nrec=length(timeout0);

  if(strcmp(test,'y'))
    it2=1;
  else
    it2=nrec;
  end

  for it=1:nspool:it2;
    timeout=timeout0(it);
    time_d=fix(timeout/86400); %days
    time_h=fix((timeout-time_d*86400)/3600);
    time_m=fix((timeout-time_d*86400-time_h*3600)/60);
    time_s=timeout-time_d*86400-time_h*3600-time_m*60;

    %out5(ivs,np) if 2D or 'S'; otherwise out5(ivs,nvrt,np)
    if(icomb==0) 
      for irank=0:nproc-1
        fname3=[base '/outputs/schout_' num2str(irank,'%04.f') '_' num2str(day) '.nc'];
        ncid3 = netcdf.open(fname3,'NC_NOWRITE');
        
        clear tmp tmp2 tmp3;
        tmp3=double(netcdf.getVar(ncid3,vid_eta,[0 it-1],[np_lcl(irank+1) 1])); 
        eta2(iplg(irank+1,1:np_lcl(irank+1)))=tmp3;

        if(ivs==1) 
          if(i23d==1) %2D
            tmp=double(netcdf.getVar(ncid3,vid5,[0 it-1],[np_lcl(irank+1) 1])); %(np_lcl)
            out5(1,iplg(irank+1,1:np_lcl(irank+1)))=tmp;
          elseif(strcmp(s_or_z,'S'))
            tmp=double(netcdf.getVar(ncid3,vid5,[lev_or_zcor-1 0 it-1],[1 np_lcl(irank+1) 1])); %(np_lcl)
            out5(1,iplg(irank+1,1:np_lcl(irank+1)))=tmp;
          elseif(strcmp(s_or_z,'Z'))
            tmp=double(netcdf.getVar(ncid3,vid5,[0 0 it-1],[nvrt np_lcl(irank+1) 1])); %(nvrt,np_lcl)
            out5(1,:,iplg(irank+1,1:np_lcl(irank+1)))=tmp;
            tmp2=double(netcdf.getVar(ncid3,vzcor,[0 0 it-1],[nvrt np_lcl(irank+1) 1])); %(nvrt,np_lcl)
            zcor(:,iplg(irank+1,1:np_lcl(irank+1)))=tmp2;
          else
            error('Unknown s_or_z');
          end
        elseif(ivs==2)
          if(i23d==1) %2D
            tmp=double(netcdf.getVar(ncid3,vid5,[0 0 it-1],[2 np_lcl(irank+1) 1])); %(ivs,np_lcl)
            out5(1:2,iplg(irank+1,1:np_lcl(irank+1)))=tmp;
          elseif(strcmp(s_or_z,'S'))
            tmp=double(netcdf.getVar(ncid3,vid5,[0 lev_or_zcor-1 0 it-1],[2 1 np_lcl(irank+1) 1])); %(ivs,np_lcl)
            out5(1:2,iplg(irank+1,1:np_lcl(irank+1)))=tmp;
          elseif(strcmp(s_or_z,'Z'))
            tmp=double(netcdf.getVar(ncid3,vid5,[0 0 0 it-1],[2 nvrt np_lcl(irank+1) 1])); %(ivs,nvrt,np_lcl)
            out5(1:2,:,iplg(irank+1,1:np_lcl(irank+1)))=tmp;
            tmp2=double(netcdf.getVar(ncid3,vzcor,[0 0 it-1],[nvrt np_lcl(irank+1) 1])); %(nvrt,np_lcl)
            zcor(:,iplg(irank+1,1:np_lcl(irank+1)))=tmp2;
          else
            error('Unknown s_or_z');
          end
        else
          error('Unknown ivs');
        end %ivs

        netcdf.close(ncid3);
      end %for irank
    else %icomb=1 (combined)
      eta2=double(netcdf.getVar(ncid0,vid_eta,[0 it-1],[np 1]));
      if(ivs==1) 
        if(i23d==1) %2D
          out5(1,:)=double(netcdf.getVar(ncid0,vid5,[0 it-1],[np 1])); %(np)
        elseif(strcmp(s_or_z,'S'))
          out5(1,:)=double(netcdf.getVar(ncid0,vid5,[lev_or_zcor-1 0 it-1],[1 np 1])); %(np)
        elseif(strcmp(s_or_z,'Z'))
          out5(1,:,:)=double(netcdf.getVar(ncid0,vid5,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
          zcor=double(netcdf.getVar(ncid0,vzcor,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
        else
          error('Unknown s_or_z');
        end
      elseif(ivs==2)
        if(i23d==1) %2D
          out5=double(netcdf.getVar(ncid0,vid5,[0 0 it-1],[2 np 1])); %(ivs,np)
        elseif(strcmp(s_or_z,'S'))
          out5=double(netcdf.getVar(ncid0,vid5,[0 lev_or_zcor-1 0 it-1],[2 1 np 1])); %(ivs,np)
        elseif(strcmp(s_or_z,'Z'))
          out5=double(netcdf.getVar(ncid0,vid5,[0 0 0 it-1],[2 nvrt np 1])); %(ivs,nvrt,np)
          zcor=double(netcdf.getVar(ncid0,vzcor,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
        else
          error('Unknown s_or_z');
        end
      else
        error('Unknown ivs');
      end %ivs
    end %icomb

    %Construct output uout(ivs,1:np) 
    if(i23d==1 || strcmp(s_or_z,'S'))
      uout=out5;
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
      h1=patch('Faces',nm','Vertices',xy00,'FaceVertexCData',uout(1,:)','FaceColor','interp','EdgeColor','none');
      colormap(jet(40));
      % Set colormap range
      caxis([0 32]); colorbar;
    else %vector
      loc_scale_x=v2(2)*0.3+v2(1)*0.7;
      loc_scale_y=-v2(4)*0.02+v2(3)*1.02;
      x_aug=[xy00(:,1)' loc_scale_x];
      y_aug=[xy00(:,2)' loc_scale_y];
      quiver(x_aug,y_aug,[uout(1,:) 1],[uout(2,:) 0]);
      text(loc_scale_x,loc_scale_y,'1 m/s');
    end %ivs

    %axis([xmin xmax ymin ymax]);
    axis([3.2e5 3.7e5 2.5e5 3.1e5]);
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
  netcdf.close(ncid0);
end %for day=stacks2
close(vidObj);
