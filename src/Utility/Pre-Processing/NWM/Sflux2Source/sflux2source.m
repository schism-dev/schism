% clear;
%----------------user inputs------------------
% make sure wdir has the following files:
%   sflux2source.prop (use xmgredit to make a ele-based prop for elements that need sources from sflux),
%   sflux/ 
%   hgrid.ll (for interpolation) 
%   hgrid.utm.26918 (for converting precip. to flux);

%example wdir='/sciclone/schism10/feiye/work/Gulf_Stream/RUN19x/Sflux2source/';
wdir='./';
sflux_files = dir([wdir '/sflux/sflux*prc_1*.nc'])
start_time_run=datenum('2017-8-4');
min_val=-0.001; max_val=100;  %values outside this range will be warned
%----------------end user inputs------------------

%---------------outputs--------------
% Please check the figure generated along with the output files 
% and see if anything is obviously wrong 
%
% Output files:
%  source_sink.in.2
%  vsource.th.2
%  msource.th.2
%------------------------------------

rho_water = 1000; %kg/m^3;

%read hgrid in meters
hgrid_name=('hgrid.utm.26918');
if exist([wdir 'hgrid.utm.mat'],'file')
    load([wdir 'hgrid.utm.mat']);%plot grid boundary
else
    display(['reading ' hgrid_name ', this can take several minutes for the first time']);
    [ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(wdir,hgrid_name,hgrid_name,0);
    save([wdir 'hgrid.utm.mat'],'ne', 'np', 'node', 'ele', 'i34', 'bndnode','open_bnds','land_bnds','ilb_island');
end

area_e=nan(ne,1); precip2flux=area_e;
for i=1:ne
    xnd_e=node(ele(i,1:i34(i)),2);
    ynd_e=node(ele(i,1:i34(i)),3);
    area_e(i)=polyarea(xnd_e(1:i34(i)),ynd_e(1:i34(i)));
    precip2flux(i) = 1/rho_water*area_e(i); % convert from  1 kg/m^2/s to m^3/s
end

%read hgrid and save a copy
hgrid_name=('hgrid.ll');
if exist([wdir 'hgrid.ll.mat'],'file')
    load([wdir 'hgrid.ll.mat']);%plot grid boundary
    figure; 
    for i=1:length(open_bnds)
        line(node(open_bnds{i},2),node(open_bnds{i},3),'LineWidth',1,'Color','k'); hold on;
    end
    for i=1:length(land_bnds)
        line(node(land_bnds{i},2),node(land_bnds{i},3),'LineWidth',1,'Color','k'); hold on;
        if ilb_island(i)==1 %island
            nd1=land_bnds{i}(end); nd2=land_bnds{i}(1);
            line(node([nd1 nd2],2),node([nd1 nd2],3),'LineWidth',1,'Color','k'); hold on;
        end
    end
else
    display(['reading ' hgrid_name ', this can take several minutes for the first time']);
    [ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(wdir,hgrid_name,hgrid_name,2);
    save([wdir 'hgrid.ll.mat'],'ne', 'np', 'node', 'ele', 'i34', 'bndnode','open_bnds','land_bnds','ilb_island');
end

lon_e=nan(ne,1); lat_e=nan(ne,1);
for i=1:ne
    lon_e(i)=mean(node(ele(i,1:i34(i)),2));
    lat_e(i)=mean(node(ele(i,1:i34(i)),3));
end

%read elements that need source
iSS=load([wdir 'sflux2source.prop']); iSS=find(iSS(:,2)==1);




nf = length(sflux_files);

fname = sflux_files(1).name;

%read basic info from the first stack of sflux
ncid = netcdf.open([wdir '/sflux/' fname],'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'time');
attval = netcdf.getAtt(ncid,varid,'base_date');
netcdf.close(ncid);
sflux_base_time = datenum([num2str(attval(1)) '-' num2str(attval(2)) '-' num2str(attval(3)) ' ' num2str(attval(4)) '0:00:00']);

lon = ncread([wdir '/sflux/' fname],'lon');
lat = ncread([wdir '/sflux/' fname],'lat');
prate = ncread([wdir '/sflux/' fname],'prate'); %kg/m^2/s
time = ncread([wdir '/sflux/' fname],'time');
dt=time(2)-time(1); nt=length(time);
dt=dt*86400; %convert to sec


fname_out=([wdir 'vsource.th.2']);
if exist(fname_out,'file')
    delete(fname_out);
end

time_stamp=0;
this_sflux_time=sflux_base_time;
for i=1:nf
    fname = sflux_files(i).name
    prate = ncread([wdir '/sflux/' fname],'prate'); %kg/m^2/s
    if min(prate(:))<-0.001
        display('Negative Precip');
    end
    for j=1:nt
        this_sflux_time=this_sflux_time+dt/86400;
        if (this_sflux_time >= start_time_run) 
            prate_interp = griddata(double(lat),double(lon),double(prate(:,:,j)),lat_e(iSS),lon_e(iSS));
            if ~(max(prate_interp(:)) < max_val & min(prate_interp(:)) > min_val) 
              display(['warning: value out of range']);
            end
            dlmwrite(fname_out,[time_stamp;max(0,prate_interp).*precip2flux(iSS)]','precision',10,'delimiter',' ','-append');
            time_stamp=time_stamp+dt;
            
            if this_sflux_time==start_time_run %diagnostic plot
                figure;
                subplot(1,3,3); 
                scatter(lon_e(iSS),lat_e(iSS),5,max(0,prate_interp),'filled'); hold on; colormap jet;colorbar;title('vsource');
            
                prate_interp = griddata(double(lat),double(lon),double(prate(:,:,j)),double(lat_e),double(lon_e));
                subplot(1,3,1); 
                contour(lon,lat,prate(:,:,j)); hold on; colormap jet; colorbar; title('sflux');              
                subplot(1,3,2); 
                [prate_sorted,II]=sort(prate_interp.*precip2flux,'descend');
                scatter(lon_e(II),lat_e(II),5,max(0,prate_interp(II)),'filled'); hold on; colormap jet;colorbar;title('vsource');
                clearvars prate_interp;
            end
        end
    end
end

fid=fopen([wdir '/source_sink.in.2'],'wt');
fprintf(fid,'%d\n',length(iSS)); %number of sources
for i=1:length(iSS)
    fprintf(fid,'%d\n',iSS(i));
end
fprintf(fid,'\n');
fprintf(fid,'%d\n',0); %number of sinks
fclose(fid);

total_nt=time_stamp/dt;
msource=zeros(total_nt,length(iSS)*2+1);
msource(:,1)=0:dt:(time_stamp-dt);
msource(:,2:length(iSS)+1)=-9999;
dlmwrite([wdir 'msource.th.2'],msource,'precision',10,'delimiter',' ');
            
    
