function []=SCHISM_TRANSECT2(base,ivs,transect_bp,varname,stacks,nspool,test)
% Authors: Joseph Zhang
% Date: Nov 2021
% Matlab function to visualize vertical transects for 3D vars (node based only) from scribed outputs.
% If there is only 1 point in the transect input, vertical 'sample' is visualized
% If there are more than 1 point in the transect input, transect animaton is produced,
% and in this case, the magnitude of vectors is plotted.
% Use inverse distance for horizontal interpolation.
% Works with mixed grid and LSC2;
% used thin layers to mimic shaved cells.
% May need to adjust some parameters inside (e.g. caxis) to get right appearance of images
% (search for 'Adjust')
% Requires get_global_info2.m (in this dir), out2d*nc, and zCoordinates*.nc
% SCHISM_TRANSECT(base,ivs,transect_bp,varname,stacks,nspool,test)
% Inputs: 
%         base: base directory (schout*.nc are in base/outputs/ and transect_bp is in base/)
%         ivs: 1 for scalar; 2 for vector
%         transect_bp: build point file name defining the transect (same projection as hgrid.gr3; depths not used)
%         varname: 'salinity' etc (node based only)
%         stacks: array of stack numbers (e.g. [2 4 5]) in the output 
%                 file names (related to time)
%         nspool: sub-sampling frequency within each stack (e.g. 1 - include all)
%         test: 'y' (only plot out 1st frame for test); 'n' (plot all frames)
% Outputs: image or transect.avi (for transect only); debug: parents.out 

close all; 
scrsz = get(0,'ScreenSize'); %screen size
figure('Position',[1 scrsz(4)*0.2 scrsz(3)/2 scrsz(4)*0.7]);

% Read transect
fid=fopen([base '/' transect_bp],'r');
c1=textscan(fid,'%d%f%f%f','headerLines',2);
%c1{1:3}: ID, x,y
fclose(fid);
xbp=c1{2}(:);
ybp=c1{3}(:);
nbp=length(xbp);
trLen=cumsum(sqrt((xbp(2:end)-xbp(1:end-1)).^2+(ybp(2:end)-ybp(1:end-1)).^2));
trLen = [0; trLen]; %trLen(1:nbp) is x-coord. of the transect

delete('transect.avi');
vidObj = VideoWriter('transect.avi');
vidObj.FrameRate = 30;  % Default 30; smaller ->slower
%vidObj.Quality = 50;    % Default 75
open(vidObj);

if(strcmp(test,'y'))
  stacks2=stacks(1); 
else
  stacks2=stacks; 
end

%axis range
xmin=min(trLen);
xmax=max(trLen)+1;
ymax=inf; 

%Get basic info
[file_nm,i23d,ne,np,nvrt,nm,xy00,vdry_e,h0,dp,kbp,vid_eta]=get_global_info2(base,varname,stacks);

% Sanity checks
if(i23d==1)
  error('2D var has no transect');
end
if(ivs==1 && nbp==1 && strcmp(test,'y'))
  error('No test sample plot for scalar; please set test to n and retry');
end

%Prep output array
out5=NaN(ivs,nvrt,np);
uout=NaN(ivs,nvrt,nbp);
zout=NaN(nvrt,nbp);

%Find parents (not Delauney tri method b/cos we have to use the weights multiple times) 
%Also avoid out of domain issue
for j=1:nbp
  iparen(j).ie=-1; %elem #
  iparen(j).i34=0; %elem type
  iparen(j).weit(1:4)=0; %weights
end %j

for ie=1:ne
  i34=4-sum(isnan(nm(:,ie))); 
  if(i34~=3 && i34 ~=4)
    i34
    error('Elem type unknown');
  end
  %in2: array of 1 (in) or 0 (out)
  in2=inpolygon(xbp,ybp,xy00(nm(1:i34,ie),1),xy00(nm(1:i34,ie),2));
  if(sum(in2)>0)
    for j=1:nbp
      if(in2(j)>0)
        iparen(j).ie=ie;
        iparen(j).i34=i34;
        %Barymetric
        clear dist idist;
        dist=sqrt((xbp(j)-xy00(nm(1:i34,ie),1)).^2+(ybp(j)-xy00(nm(1:i34,ie),2)).^2)+1.e-20; %avoid /0
        idist=1./dist;
        sum2=sum(idist(1:i34));
        iparen(j).weit(1:i34)=idist(1:i34)/sum2;
      end %if
    end %for j

    %Check if done
    for j=1:nbp
      itmp5(j)=iparen(j).ie;
    end
    if(length(find(itmp5(1:nbp)<=0))==0) 
      disp('finished searching for parents at elem #:');
      ie
      break; 
    end
  end %if
end %for ie

%Output
fid2=fopen('parents.out','w');
ifl=0;
for ii=1:nbp
  fprintf(fid2,'Point #%d: %d %d %f %f %f %f %f\n',ii,iparen(ii).ie,iparen(ii).i34,...
          iparen(ii).weit(1:4),sum(iparen(ii).weit(1:4)));
  if(iparen(ii).ie<0)
    disp(['Parent not found for point #:' num2str(ii)]);
    ifl=1;
  end
  ie=iparen(ii).ie;
end %ii
if(ifl==1); error('Check parents.out for errors'); end;
fclose(fid2);

count=0;
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

  dtout=timeout0(2)-timeout0(1);
  time_min=(stacks2(1)-1)*nrec*dtout/86400; %for axis (days)
  time_max=stacks2(end)*nrec*dtout/86400;

  if(strcmp(test,'y'))
    it2=1;
  else
    it2=nrec;
  end

  %Open out2d and zcor
  fname_2d=[base,'/outputs/out2d_' num2str(day) '.nc'];
  ncid_2d = netcdf.open(fname_2d,'NC_NOWRITE');

  fname_z=[base,'/outputs/zCoordinates_' num2str(day) '.nc'];
  ncid_z = netcdf.open(fname_z,'NC_NOWRITE');
  vid_z=netcdf.inqVarID(ncid_z,'zCoordinates');

  for it=1:nspool:it2;
    count=count+1;
    % Time info
    timeout=timeout0(it); %sec
    time_d=fix(timeout/86400); %days
    time_h=fix((timeout-time_d*86400)/3600);
    time_m=fix((timeout-time_d*86400-time_h*3600)/60);
    time_s=timeout-time_d*86400-time_h*3600-time_m*60;

    zcor=double(netcdf.getVar(ncid_z,vid_z,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
    idry_e=netcdf.getVar(ncid_2d,vdry_e,[0 it-1],[ne 1]); %(ne)

    out5(1,:,:)=double(netcdf.getVar(ncid1,vid1,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
    if(ivs==2)
      out5(2,:,:)=double(netcdf.getVar(ncid2,vid2,[0 0 it-1],[nvrt np 1])); %(nvrt,np)
    end %ivs

    %Output array: uout(ivs,nvrt,nbp) 
    for i=1:nbp
      ie=iparen(i).ie;
      %Init for below bottom etc
      uout(:,:,i)=nan;
      zout(:,i)=nan;
      if(idry_e(ie)~=0)
        kbp_e(i)=0;
        continue; 
      end
     
      %Wet elem
      i34=iparen(i).i34;
      weit(1:i34)=iparen(i).weit(1:i34);
      nds(1:i34)=nm(1:i34,ie);
      kbp_e(i)=min(kbp(nds(1:i34)));
      for k=kbp_e(i):nvrt
        zout(k,i)=0;
        uout(:,k,i)=0;
        for m=1:i34
          zout(k,i)=zout(k,i)+zcor(max([k kbp(nds(m))]),nds(m))*weit(m);
          uout(:,k,i)=uout(:,k,i)+out5(1:ivs,max([k kbp(nds(m))]),nds(m))*weit(m);
        end %for m
      end %for k
    end %for i=1:nbp

    %Pad extra thin layers for pcolor
    kbout=min(kbp_e(1:nbp));
    if(kbout<1); kbout=1; end;
%    kbout
%    kbp_e(1:nbp)
    for i=1:nbp
      for j=1:ivs
        uout(j,kbout:kbp_e(i)-1,i)=uout(j,max([1 kbp_e(i)]),i);
      end %j
      clear tmp;
      tmp=(kbout-kbp_e(i)):-1;
      zout(kbout:kbp_e(i)-1,i)=zout(max([1 kbp_e(i)]),i)+tmp*1.e-7;
    end %for i=1:nbp

    %Define ymin
    if(count==1)
      ymin=min(min(zout))-0.5;
    end

    % plot
    if(nbp==1) %sample
      %Contruct matrix for scalars (plot later; inside it loop to conserve memory)
      sout(:,count)=uout(1,:,1);
      if(ivs==2); tout(:,count)=uout(2,:,1); end;
      sxout(:,count)=timeout/86400*ones(nvrt,1); %x coord.
      syout(:,count)=zout(:,1); %y coord.
%      continue;
    else %transect
      hold on;
      axis([xmin xmax ymin ymax]);
      v2=axis;
      % Write time stamp info
      loc_info_x=(v2(2)+v2(1))/2;
      loc_info_y=v2(4)*0.95+v2(3)*0.05;
%      text(loc_info_x,loc_info_y,{'Time (DD:HH:MM:SS)'; num2str([time_d time_h time_m time_s])});
      title(['Time (DD:HH:MM:SS): ' num2str([time_d time_h time_m time_s])]);

      vmag=squeeze(uout(1,:,:));
      if(ivs==2)
        vmag=squeeze(sqrt(uout(1,:,:).^2+uout(2,:,:).^2));
      end
      h=pcolor(trLen,zout,vmag);
      set(h,'EdgeColor','none','FaceColor','interp');
      set(gcf,'Color',[1 1 1]);
      colormap(jet(40));
      %Adjust
      caxis([0 30]); colorbar;
      xlabel('Along transect distance (m)');
      ylabel('z (m)');

      %Plot out transect grid here if desired
      %plot(trLen,zout','k.');

      set(gca,'nextplot','replacechildren');
      currFrame = getframe(gcf);
      writeVideo(vidObj,currFrame);
      if(day ~= stacks2(end) || it ~=it2)
        clf; %clear figure to avoid overlay
      end
    end %nbp
  end %it

  netcdf.close(ncid1);
  netcdf.close(ncid_2d);
  netcdf.close(ncid_z);
  if(ivs==2); netcdf.close(ncid2); end;
end %for day=stacks2

close(vidObj);

%Finish plotting for sample
%sout(nvrt,count) (count is # of time steps); tout (for vectors)
if(nbp==1)
  axis([time_min*0.99 time_max*1.01 ymin ymax]);
  hold on;

  if(ivs==1) 
    h=pcolor(sxout,syout,sout);
    set(h,'EdgeColor','none','FaceColor','interp');
    colormap(jet(40));
    %Adjust
    caxis([0 33]); colorbar;
  else %vector
    quiver([sxout; sxout(end,:)],[syout; (ymin*0.98+ymax*0.02)*ones(1,count)], ...
    [sout; 1*ones(1,count)],[tout; 0*ones(1,count)],0,'k');
    text(time_min*0.7+time_max*0.3,ymin*0.98+ymax*0.02,'1 m/s');
  end %ivs
  set(gcf,'Color',[1 1 1]);
  xlabel('Time (days)');
  ylabel('z (m)');
end %sample plot
