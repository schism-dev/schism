%Plot output from gen_vqs.f90 (vgrid_master.out, transect*.out)
close all; clear all;
z_m=load('vgrid_master.out'); %m,nv(m),hsm,z_mas()
z1=load('transect1.out'); %i,kbp,x,y,transect_dis,dp,z-coor

np_m=size(z_m,1)
nv_m=size(z_m,2)-3
zcor_m=z_m(:,4:end);
for i=1:np_m
  kbp_m=z_m(i,2);
  zcor_m(i,kbp_m+1:end)=nan;
end %for

np=size(z1,1)
nvrt=size(z1,2)-6
zcor1=z1(:,7:end);
for i=1:np
  kbp=z1(i,2);
%  zcor1(i,kbp+1:end)=nan;
end %for

subplot(2,1,1); hold on;
plot(z_m(:,1),zcor_m,'k-',z_m(:,1),-z_m(:,3),'r.');
for i=1:np_m
  plot(z_m(i,1)*ones(nv_m,1),zcor_m(i,:),'k');
end %for i
title('Master grid');
xlabel('Grid #');
subplot(2,1,2); hold on;
plot(z1(:,5),zcor1(:,:),'k-',z1(:,5),-z1(:,6),'r.');
for i=1:np
  plot(z1(i,5)*ones(nvrt,1),zcor1(i,:),'k');
end %for i
title('Transect before adjustment (transect1)');
xlabel('Along transect distance (m)');
