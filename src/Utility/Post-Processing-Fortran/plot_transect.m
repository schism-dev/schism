%Plot (scalar) transect output from read_output*_transect.f90. Can produce either transect animation or 'sample' view 
%for each station (time-depth).
%Inputs: transect.out, transect_grd.r0, transect_grd.z0. The transect grid must be fixed
%        In time. Inputs at the beginning.
%Edit plot appearance below.
clear all; close all;

iflag_transect=0; %0: sample for each station; 1: transect animation
run='outputs/'; %path to data
var=load([run '/transect.out']); %time(sec),out3(1:nxy*ntranz,it) [vertical direction first]
r0=load([run '/transect_grd.r0']); %along transect distance
z0=load([run '/transect_grd.z0']); % z (m)
start_year=2019;
start_mon=9;
start_day=1;

nr=length(r0); nz=length(z0);
nrz=nr*nz;
for i=1:nr
  X(i,1:nz)=r0(i);
  for j=1:nz
    Z(i,j)=z0(j);
  end %for j 
end %for i

ntime=size(var,1)/nrz;

if(iflag_transect==0) %sample
  timeout=nan(ntime,nz);
  zout=timeout;
  var3=nan(ntime,nz,nr);;
  for it=1:ntime
    timeout(it,1:nz)=var(it*nrz,1);
    zout(it,1:nz)=z0;
    var3(it,1:nz,1:nr)=reshape(var((it-1)*nrz+1:it*nrz,2),nz,nr);
  end %it
  var3(find(var3<-900))=nan;

  ifig=0;
  for ista=1:nr
    ifig=ifig+1;
    figure(ifig);
    pcolor(timeout/86400,zout,squeeze(var3(:,:,ista)));
    title(['Station #' num2str(ista)]);
    %ylim([-40 0]); %set individual y axis here
    xlabel('Days');
    ylabel('Z (m)');
    shading interp;
    caxis([0 20]); 
    colormap(jet(40));
    colorbar;
  end %ista
else %transect
  avi_out = avifile('anim_transect_T.avi','FPS',5); %output movie to slab.avi
  figure(1);
  for it=1:ntime
    time=var(it*nrz,1);
    date=datestr(datenum(start_year,start_mon,start_day)+time/86400,31);
    var2=var((it-1)*nrz+1:it*nrz,2);
    var2(find(var2<-900))=nan;
    var3=reshape(var2,nz,nr);
  
    subplot(2,1,1);
    pcolor(X,Z,var3');
    shading interp;
    ylim([-250 0]);
    caxis([-2 30]);
    colormap(jet(40));
    colorbar;
    title([date ' ; T']);
    % Add image to avi file
    frame = getframe(gcf);
    avi_out=addframe(avi_out,frame);
    clf;
      
    clear var2 var3;
  end %for it
  avi_out=close(avi_out);
end %sample or transect
