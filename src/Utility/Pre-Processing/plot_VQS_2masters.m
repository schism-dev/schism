%Plot output from gen_vqs_2masters.f90 (vgrid_master*.out, transect1.out)
close all; clear all;
z1=load('transect1.out'); %i,kbp,x,y,transect_dis,dp,set_flag,z-coor
n_master=2; %# of masters

%Masters
figure(1);
for m=1:n_master
  clear z_m zcor_m;
  z_m=load(['vgrid_master' num2str(m) '.out']); %m,nv(m),hsm,z_mas()
  np_m=size(z_m,1)
  nv_m=size(z_m,2)-3
  zcor_m=z_m(:,4:end);
  for i=1:np_m
    kbp_m=z_m(i,2);
    zcor_m(i,kbp_m+1:end)=nan;
  end %for i

  subplot(2,1,m); hold on;
  plot(z_m(:,1),zcor_m,'k-',z_m(:,1),-z_m(:,3),'r.');
  for i=1:np_m
    plot(z_m(i,1)*ones(nv_m,1),zcor_m(i,:),'k');
  end %for i
  title(['Master grid #' num2str(m)]);
  xlabel('Grid #');
end %m

%Transect plot
figure(2);
np=size(z1,1)
nvrt=size(z1,2)-7
zcor1=z1(:,8:end);
%for i=1:np
%  kbp=z1(i,2);
%  zcor1(i,kbp+1:end)=nan;
%end %for
plot(z1(1,5),zcor1(1,:),'k-',z1(1,5),-z1(1,6),'r.'); %in case 1 pt is 'isolated'
start=1; %Plot each segment of build pts
for i=1:np
  if(i>1); if(abs(z1(i,7)-z1(i-1,7))>1.e-3);
    plot(z1(start:i-1,5),zcor1(start:i-1,:),'k-',z1(start:i-1,5),-z1(start:i-1,6),'r.');
    start=i; %reset starting index
  end; end;

  plot(z1(i,5)*ones(nvrt,1),zcor1(i,:),'k');
end %for i
%Plot last seg
plot(z1(start:np,5),zcor1(start:np,:),'k-',z1(start:np,5),-z1(start:np,6),'r.');

title('Transect before adjustment (transect1)');
xlabel('Along transect distance (m)');
