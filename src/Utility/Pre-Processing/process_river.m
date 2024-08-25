%Extract river flows and prep vsource.th etc (up to 2018) from (Dai & Trenberth). Needs read_gr3.m
%Script will combine duplicated source elem's. River mouth locations are used in bp.
%Inputs: out_river (gredit extract from coastal-stns-byVol-updated-oct2007.bp, commented out 1st 2 lines. 
%        lon in [-180,180]); 
%        hgrid.gr3 (lon in [-180,180])
%        coastal-stns-Vol-monthly.updated-May2019.nc (Dai&Trenberth via ucar)
%        dates below
%Outputs: source_sink.in, vsource.th, msource.th (-9999 for T; S specified by user). Time step is 1 day and may need pad 
%         extra record for last model step

clear all; close all;

%Output info
river_S=30;
start_yr=2017; end_yr=2018;
start_mon=12; end_mon=12;
start_day=1; end_day=31; %GMT

%Location
riv=load('out_river'); %id,lon,lat,rank
nriv_out=size(riv,1);

ncid0 = netcdf.open('coastal-stns-Vol-monthly.updated-May2019.nc','NC_NOWRITE');
vid=netcdf.inqVarID(ncid0,'time'); 
time= double(netcdf.getVar(ncid0, vid)); %YYYYMM; (1:nmonths)
vid=netcdf.inqVarID(ncid0,'FLOW'); 
monthly_flow= double(netcdf.getVar(ncid0, vid)); %m^3/s; (nrivers,nmonths)
[nrivers nmonths]=size(monthly_flow);
%starting yr, month
nc_start_yr=fix(time(1)/100);
nc_start_mon=rem(time(1),100);

%Calc monthly mean to fill in junks
icount(1:nrivers,1:12)=0;
mean_flow(1:nrivers,1:12)=0;
for i=1:nmonths
  yr=fix(time(i)/100);
  mon=fix(time(i)-100*yr);
  for j=1:nrivers
    if(monthly_flow(j,i)>0)
      mean_flow(j,mon)=mean_flow(j,mon)+monthly_flow(j,i);
      icount(j,mon)=icount(j,mon)+1;
    end
  end %for j
end %for i
mean_flow=mean_flow./icount;

%Save time series at all rivers first for all output months. Replace junks with clima
nmon=(end_yr-start_yr)*12+end_mon-start_mon+1; %# of output months
timeout(1:nmon)=0; %time pts in days
out_flow(1:nrivers,1:nmon)=0;
yr_now=start_yr;
mon_now=start_mon;
day_now=start_day;
for i=1:nmon
  timeout(i)=datenum(yr_now,mon_now,1)-datenum(start_yr,start_mon,start_day); %days
  record=(yr_now-nc_start_yr)*12+mon_now-nc_start_mon+1; %time record # in .nc
  for j=1:nrivers
    if(monthly_flow(j,record)>0)
      out_flow(j,i)=monthly_flow(j,record); %flow
    else
      out_flow(j,i)=mean_flow(j,mon_now); %flow
    end
  end %for j

  mon_now=mon_now+1;
  if(mon_now>12) 
    mon_now=mon_now-12; 
    yr_now=yr_now+1;
  end
end %for i

timeout2=0:timeout(end); %final output is daily
nstep=length(timeout2);
out_flow2(1:nriv_out,1:nstep)=0;
for j=1:nriv_out
  iriv=riv(j,4); %river rank as index into original nc file
  out_flow2(j,:)=interp1(timeout,out_flow(iriv,:),timeout2);
end %for j
out_flow2(find(out_flow2<0))=0;
out_flow2(isnan(out_flow2))=0;
%if(sum(sum(find(out_flow2<0)))>0)
%  error('Negative flow');
%end

%Read in hgrid.gr3
[xnd,ynd,dp,i34,nm]=read_gr3('hgrid.gr3');
np=length(xnd);
ne=length(i34);

%Elem center (not accurate across dateline)
for i=1:ne
  xctr(i)=sum(xnd(nm(1:i34(i),i)))/double(i34(i));
  yctr(i)=sum(ynd(nm(1:i34(i),i)))/double(i34(i));
end %for i

%Find nearest elem
for i=1:nriv_out
  disp(['doing river:' num2str(i)]);
  dist2=(xctr-riv(i,2)).^2+(yctr-riv(i,3)).^2;
  [minD(i),ie_river(i)]=min(dist2);
end %for i

%Combine redundant elem
nriv_out2=0; %final # of output rivers
elem_to_rank(1:ne)=fix(0.); %point to new rank from 1st found river 
for i=1:nriv_out
  ie=ie_river(i);
  if(elem_to_rank(ie)==0) %new river
    nriv_out2=nriv_out2+1;
    out_flow3(nriv_out2,1:nstep)=out_flow2(i,:);
    elem_to_rank(ie)=fix(nriv_out2);
    rank_to_elem(nriv_out2)=fix(ie);
  else %add to existing
    out_flow3(elem_to_rank(ie),1:nstep)=out_flow3(elem_to_rank(ie),1:nstep)+out_flow2(i,:);
  end
end %for i

%Output
out_char='%f';
out_char2='%f'; %msource.th
for i=1:nriv_out2
  out_char=[out_char ' %f'];
  out_char2=[out_char2 ' %f %f'];
end

fid=fopen('vsource.th','w');
fid2=fopen('msource.th','w');
fprintf(fid,[out_char '\n'],[timeout2*86400; out_flow3]);
fprintf(fid2,[out_char2 '\n'],[timeout2*86400; -9999*ones(nriv_out2,nstep); river_S*ones(nriv_out2,nstep)]);
fclose(fid); fclose(fid2);

fid=fopen('source_sink.in','w');
fprintf(fid,'%d\n',nriv_out2);
for i=1:nriv_out2
  fprintf(fid,'%d\n',rank_to_elem(i));
end
fprintf(fid,'\n %d',0);
fclose(fid);
