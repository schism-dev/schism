%Plot (scalar) transect output from read_output*_transect.f90
%Inputs: transect.out, transect_grd.r0, transect_grd.z0. The transect grid must be fixed
%        In time.
clear all; close all;
run='RUN22a';
var=load(['../' run '/transect.out']); %time(sec),out3(1:nxy*ntranz,it) [vertical direction first]
r0=load(['../' run '/transect_grd.r0']); %along transect distance
z0=load(['../' run '/transect_grd.z0']); % z (m)
start_year=2011;
start_mon=7;
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
