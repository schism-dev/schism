function [xnd,ynd,dp,i34,nm]=read_gr3(fname)
%xnd,ynd,dp,i34,nm]=read_gr3(fname)
%e.g. read_gr3({'hgrid.gr3'},
% fname: .gr3 name
%Outputs: xnd(np),ynd,dp,i34,nm(4,ne)

  fid=fopen(fname,'r');
  char=fgetl(fid);
  tmp1=str2num(fgetl(fid));
  fclose(fid);
  
  ne=fix(tmp1(1));
  np=fix(tmp1(2));
  
  fid=fopen(fname,'r');
  c1=textscan(fid,'%d%f%f%f',np,'headerLines',2);
  fclose(fid);
  fid=fopen(fname,'r');
  c2=textscan(fid,'%d%d%d%d%d%d',ne,'headerLines',2+np);
  fclose(fid);
  
  xnd=c1{2}(:);
  ynd=c1{3}(:);
  dp=c1{4}(:);
  i34=c2{2}(:);

  %Make lon in [0,360)
%  if(isphere~=0)
%    indx=find(x<0);
%    x(indx)=x(indx)+360;
%  end
  
  nm(1:4,1:ne)=nan;
  for i=1:ne
    for j=1:i34(i)
      nm(j,i)=fix(c2{j+2}(i));
    end %for j

    %Check discontinuity across prime meridian
%    if(isphere~=0)
%      ifl=0; %flag
%      for j=1:i34(i)
%        n1=nm(i,j);
%        j2=j+1;
%        if(j==i34(i)); j2=j2-i34(i); end;
%        n2=nm(i,j2);
%        if(abs(x(n1)-x(n2))>180.)
%          ifl=1; break;
%        end
%      end %for j
%
%      if(ifl>0) %mask out this elem
%        nm(i,:)=nan;
%      end
%    end %isphere/
  end %for i
