%Read tidal constants and output tidal potential
clear all; close all;
load tide_fac_constants.mat;
nline=length(const.name);
aa=cellstr(const.name); %convert to cell string for easy write
fid=fopen('tidal_const.out','w');
for i=1:nline
  %freq converted to rad/sec
  fprintf(fid,'%4s %.10f %f %f %f\n',aa{i},const.freq(i)*2*pi/3600,const.doodsonamp(i), ...
          const.doodsonspecies(i),const.potentialamp(i));
end %for
fclose(fid);
