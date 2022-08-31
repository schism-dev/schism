function []=sflux2source(varargin)
%----------------user inputs------------------
% make sure wdir has the following files:
%   sflux2source.prop
%     (use xmgredit to make an ele-based prop file for elements (prop value=1) that need sources from sflux,
%     if you don't have xmgredit, use the alternative method shown in /Sflux2Source/Sample_Prop/)
%   sflux/
%   hgrid.ll (for interpolation)
%   hgrid.utm.26918 (for converting precip. to flux);
%
% Important! To speed up the interpolation, the script assumes
%   the lon/lat does not change among "sflux_fliles".
%   Different interpolation method will return slightly different vsources.
%   e.g., the inverse distance weighted method generally produces worse patterns than linear interp

if nargin==0
    %Sample inputs. You can change the variables here if you want to run this script directly.
    %1) it is recommended to copy (cp -rL) the entire script folder (Sflux2Source) into your run dir to keep a record of the scripts used
    wdir='./';
    %2) look for sflux files using the 'dir' function. Depending on your needs,
    %   you may change the path, e.g., sflux_files=dir('../sflux/sflux*prc_1*.nc'),
    %   or set the sflux dataset number, e.g., sflux_files = dir([some_dir sflux/sflux*prc_2*.nc']);
    sflux_files = dir([wdir '/sflux/sflux*prc_1*.nc']);
    %3) start time of the SCHISM run
    start_time_run=datenum('2018-6-17 00:00:00');
    %4) values outside this range will be warned
    min_val=-0.001; max_val=100;  
    %5) interp method: 
    i_interp=1;  % 0: fastest, 8-point inverse distance weighted, assuming lon/lat does not change in sflux files
                 % 1 (default): linear interp with scatteredInterpolant, assuming lon/lat does not change in sflux files, 4 times slower than 0
                 % 2: linear interp with shape functions, assuming lon/lat does not change in sflux files, x times slower than 0; needs more test
elseif nargin==6
    %Run the script like a function
    wdir=varargin{1,1};
    sflux_files=varargin{1,2};
    start_time_run=varargin{1,3};
    min_val=varargin{1,4};
    max_val=varargin{1,5};
    i_interp=varargin{1,6};
else
    disp('wrong number of input arguments');
end
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

disp(['converting sflux files from <' sflux_files(1).name '> to <' sflux_files(end).name '>'])

rho_water = 1000; %kg/m^3;

%read hgrid in meters
hgrid_name=('hgrid.utm.26918');
display(['reading ' hgrid_name]);
[ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(wdir,hgrid_name,hgrid_name,0);
display(['done reading ' hgrid_name]);
   
display(['calculating element area based on ' hgrid_name]);
area_e=nan(ne,1); precip2flux=area_e;
for k=3:4
    nn = find(i34==k);
    area_x=reshape(node(ele(nn,1:k)',2),k,length(nn));
    area_y=reshape(node(ele(nn,1:k)',3),k,length(nn));
    area_e(nn)=polyarea(area_x,area_y);
end
precip2flux = 1/rho_water*area_e; % convert from  1 kg/m^2/s to m^3/s

%read hgrid and save a copy
hgrid_name=('hgrid.ll');
display(['reading ' hgrid_name]);
[ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(wdir,hgrid_name,hgrid_name,0);
display(['done reading ' hgrid_name]);


display(['calculating element center based on ' hgrid_name]);
node=[node; [nan nan nan nan]]; % add a dummy node at the end
ele(isnan(ele))=np+1;
lon_e=nanmean(reshape(node(ele',2),4,ne)',2);
lat_e=nanmean(reshape(node(ele',3),4,ne)',2);
node=node(1:end-1,:); ele(ele==np+1)=nan; %restore node and ele

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
dt=round(dt*86400); %convert to sec

%calculate weights for interpolation
if i_interp==0 % assuming lon/lat of sflux*.nc does not change
    disp('calculating weights for spatial interpolation')
    n_inv_points = 4
    X = [double(lat(:)),double(lon(:))];
    Mdl = createns(X,'Distance','euclidean');
    Q = [lat_e(iSS),lon_e(iSS)];
    IdxNN = knnsearch(Mdl,Q,'K',n_inv_points);
    dist = nan*IdxNN;
    for k=1:n_inv_points
        dist(:,k)=( (Q(:,1)-X(IdxNN(:,k),1)).^2 + (Q(:,2)-X(IdxNN(:,k),2)).^2 ).^0.5;
    end
    inv_weight = 1./dist ./ sum(1./dist, 2);
elseif i_interp==1
    this_prate=prate(:,:,1);
    F_interp=[]; seq3=[1 2 3 1 2];
    F_interp = scatteredInterpolant(double(lat(:)),double(lon(:)),double(this_prate(:)));                
elseif i_interp==2
    disp('calculating weights for linear interpolation')
    this_prate=prate(:,:,1);
    t=delaunayn([double(lat(:)),double(lon(:))]);
    f=double(this_prate(:));
    i_ele = tsearchn([double(lat(:)),double(lon(:))],t,[lat_e(iSS),lon_e(iSS)]);
    i_area_cor=-ones(length(i_ele),3);

    for k=1:3
        area_lat0 = [lat_e(iSS) lat(t(i_ele,seq3(k+1:k+2)))]';
        area_lon0 = [lon_e(iSS) lon(t(i_ele,seq3(k+1:k+2)))]';
        area_lat = lat(t(i_ele,:))'; area_lon = lon(t(i_ele,:))';
        i_area_cor(:,k) = abs(polyarea(area_lat0,area_lon0)./polyarea(area_lat,area_lon));
    end

    II=find(i_area_cor<0 | i_area_cor>1);
    if ~isempty(II)
        disp('error in area coordinates')
    end
    disp('done calculating weights for linear interpolation')
end


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
        disp('Negative Precip');
    end
    for j=1:nt
        this_sflux_time=this_sflux_time+double(dt/86400);
        if (this_sflux_time >= start_time_run)
           
            this_prate=prate(:,:,j);
            if i_interp==0
                prate_1d = this_prate(:);
                prate_interp = sum(prate_1d(IdxNN) .* inv_weight, 2);
            elseif i_interp==1
                F_interp.Values=double(this_prate(:));
                prate_interp = F_interp(lat_e(iSS),lon_e(iSS));
            elseif i_interp==2
                tmp=this_prate(:);
                this_prate_inele = tmp(t(i_ele,:));
                prate_interp = sum(this_prate_inele .* i_area_cor, 2);
            end
           
            if ~(max(prate_interp(:)) < max_val && min(prate_interp(:)) > min_val)
              disp(['warning: value out of range']);
              max(prate_interp(:))
              min(prate_interp(:))
            end
            dlmwrite(fname_out,[time_stamp;max(0,prate_interp).*precip2flux(iSS)]','precision',10,'delimiter',' ','-append');
            time_stamp=time_stamp+dt;
                        
            if this_sflux_time==start_time_run %diagnostic plot
                figure;
                subplot(1,3,3); 
                scatter(lon_e(iSS),lat_e(iSS),5,max(0,prate_interp),'filled'); hold on; colormap jet;colorbar;title('vsource, selected region');
                
                subplot(1,3,1); 
                contour(lon,lat,prate(:,:,j)); hold on; colormap jet; colorbar; title('sflux'); 
                
                subplot(1,3,2); 
                if i_interp==0
                    Q0 = [lat_e,lon_e];
                    IdxNN0 = knnsearch(Mdl,Q0,'K',n_inv_points);
                    dist0 = nan*IdxNN0;
                    for k=1:n_inv_points
                        dist0(:,k)=( (Q0(:,1)-X(IdxNN0(:,k),1)).^2 + (Q0(:,2)-X(IdxNN0(:,k),2)).^2 ).^0.5;
                    end
                    inv_weight0 = 1./dist0 ./ sum(1./dist0, 2);
                    this_prate0=prate(:,:,j);
                    prate_1d0 = this_prate0(:);
                    prate_interp = sum(prate_1d0(IdxNN0) .* inv_weight0, 2);                 
                    scatter(lon_e,lat_e,5,max(0,prate_interp),'filled'); hold on; colormap jet;colorbar;title('vsource, whole domain');
                elseif i_interp==1
                    this_prate=prate(:,:,j);
                    F_interp.Values = double(this_prate(:));
                    prate_interp = F_interp(double(lat_e),double(lon_e));
                    scatter(lon_e,lat_e,5,max(0,prate_interp),'filled'); hold on; colormap jet;colorbar;title('vsource, whole domain');
                elseif i_interp==2
                    % do nothing, skip subplot 2
                end

                % [prate_sorted,II]=sort(prate_interp0.*precip2flux,'descend');
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

msource=zeros(2,length(iSS)*2+1);
msource(:,1)=[0 time_stamp*1.1];
msource(:,2:length(iSS)+1)=-9999;
dlmwrite([wdir 'msource.th.2'],msource,'precision',15,'delimiter',' ');
