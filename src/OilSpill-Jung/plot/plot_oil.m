clear all; close all;
pathdata = partpath('../run/spill.out');
%pathdata is a struc:
%   nsteps: # of time steps;
%   npar: # of particles;
%   pathdata.time(1:nsteps): time (sec)
%   pathdata.x(1:npar,1:nsteps): x coordinates
%   pathdata.y(1:npar,1:nsteps): y coordinates
%   pathdata.z(1:npar,1:nsteps): z coordinates (relative to MSL after moving; intial value before moving)
%   pathdata.el(1:nsteps): local elevtion at 1st particle;
%   pathdata.wx(1:nsteps): local wind velocity in x at 1st particle;
%   pathdata.wy(1:nsteps): local wind velocity in y at 1st particle;
[npar,nsteps]=size(pathdata.x);
for it=1:20 %nsteps
  subplot(4,5,it);
  plot(pathdata.x(:,it),pathdata.y(:,it),'ob','MarkerFaceColor','b');
  title(['Time (days)=' num2str(pathdata.time(it)/86400)]);
end %for
