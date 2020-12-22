function []=plot_gr3(fname,caxis_min,caxis_max,num_columns,isphere)
%plot_gr3(fname,caxis_min,caxis_max,num_columns) 
%Plot depths in .gr3 (tri-quad) in matlab
% where fname is a cell array (e.g. {'a','b'}
%caxis_min,caxis_max are min/max used in caxis
%num_columns: # of columna in subplot
%isphere: shperical grid option. If /=0, assumes lon discontinuity @dateline, 
%         and will recast lon to [0,360) and mask out discontinuity @prime meridian

close all;
nfiles=length(fname);
nrows=ceil(nfiles/num_columns);
figure(1);
set(gcf,'Color',[1 1 1]);
for ifile=1:nfiles
  fid=fopen(fname{ifile},'r');
  char=fgetl(fid);
  tmp1=str2num(fgetl(fid));
  fclose(fid);
  
  ne=fix(tmp1(1));
  np=fix(tmp1(2));
  
  fid=fopen(fname{ifile},'r');
  %Change here if there are >1 'depths'
  c1=textscan(fid,'%d%f%f%f',np,'headerLines',2);
  fclose(fid);
  fid=fopen(fname{ifile},'r');
  c2=textscan(fid,'%d%d%d%d%d%d',ne,'headerLines',2+np);
  fclose(fid);
  
  x=c1{2}(:);
  y=c1{3}(:);
  bathy=c1{4}(:);
  i34=c2{2}(:);

  %Make lon in [0,360)
  if(isphere/=0)
    indx=find(x<0);
    x(indx)=x(indx)+360;
  end
  
  nm(1:ne,1:4)=nan;
  for i=1:ne
    for j=1:i34(i)
      nm(i,j)=fix(c2{j+2}(i));
    end %for j

    %Check discontinuity across prime meridian
    if(isphere/=0)
      ifl=0; %flag
      for j=1:i34(i)
        n1=nm(i,j);
        j2=j+1;
        if(j==i34(i)); j2=j2-i34(i); end;
        n2=nm(i,j2);
        if(abs(x(n1)-x(n2))>180.)
          ifl=1; break;
        end
      end %for j

      if(ifl>0) %mask out this elem
        nm(i,:)=nan;
      end
    end %isphere/
  end %for i
  
  subplot(nrows,num_columns,ifile); hold on;
  %Plot with grid on
  %patch('Faces',nm(:,1:4),'Vertices',[x y],'FaceVertexCData',bathy,'FaceColor','interp'); 
  patch('Faces',nm(:,1:4),'Vertices',[x y],'FaceVertexCData',bathy,'FaceColor','interp','EdgeColor','none');
  colormap(jet);
  caxis([caxis_min caxis_max]);
  if(ifile==nfiles); colorbar; end;
  %axis([5e5 6e5 4.13e6 4.24e6]);
end %for ifile
