%Author: Huy Tran (huy2013uq@gmail.com)
clear all; close all; 
%% CREATE bctides.in for SCHISM in one-click
% Creating bctides.in for SCHISM becomes easier (purely in MATLAB) & thanks to  OceanMesh2D's TEAM
% https://github.com/CHLNDDEV/OceanMesh2D
% bctides.in  can be created in least than a minute with ONE-CLICK as I
% have put all functions that need to create this file inside this script.
% Main functions are:
%    *  extract_open_nodes:  uses readfort14 in  OceanMesh2D github
%    ** gen_harm_FES_SCHISM: a modifed version of  gen_harm_FES.m in SCHISM
%    *** tide_fac: in  OceanMesh2D github  needs tide_fac_constants.mat
%   (using this the user does not need to run tidal_fac.f)

% INPUTS:
%            *   FES2014: downloaded and placed in a directory
%               (e.g. C:\FES_TIDE_2014)
%            ** hgrid.gr3 (lon/lat) to get the open boundary nodes
%            *** tide_fac_constants.mat (provided)
%            (https://github.com/CHLNDDEV/OceanMesh2D/blob/master/utilities/tide_fac_constants.mat)

% OUTPUTS: 
%           * bctides.in: will be created in a directory e.g. C:\FES_TIDE_2014\INPUTS

% %%%%%%%% INPUT INFORMATION %%%% %%%%%%%%%%%%
%% 1- where to save bctides.in?
savedirs='C:\FES_TIDE_2014\INPUTS\';
bounodes=extract_open_nodes(strcat(savedirs,'hgrid.gr3')); 
%% 2- Specify time to run model
t_s='01-Jan-2008 01:00:36'; % Start
t_e='01-Jan-2009 01:00:36'; % Finish
%% 3- Specify how many tidal constituents you expect to run?
 incnstit = {'M2','S2','N2','K2','K1','O1','Q1','P1'};
 
% incnstit={'T2','SSA','SA','S2','S1','R2','Q1',...
% 'P1','O1','NU2','N2','MU2','MSF',...
% 'MM','MF','M3','M2',...
% 'L2','K2','K1','J1','EPS2','2N2'};

%% 4- Select flag options (e.g run with tidal elevation only)
iflag=2;            % 1: only tidal elevations; 2: elevations + uv velocity;
ops=[3 0 0 0];   % refer to SCHSIM manuals for details.

%%%% %%%%%%%%%%%% YOU MAY NOT NEED TO EDIT  LINES BELOW %%%% %%%%%%%%%%%%
obj = tide_fac(cell2mat(bounodes(:)),t_s,t_e,incnstit);
%% Write tide_fact.out
fid = fopen(strcat(savedirs, 'bctides.in'), 'w' ) ;
fprintf( fid, '%s %s %d %s  \n', t_s,'!', datenum(t_e)-datenum(t_s),'day'); %
fprintf( fid, '%d %s \n',  obj.f15.ntif , '1000. ntip') ; %

for k = 1: obj.f15.ntif
lastr={obj.f15.tipotag.name};
xa(k)=str2double(lastr{1,k}(end));
xa(isnan(xa))=0; xa(xa>2)=0;
fprintf( fid, '%s \n',  obj.f15.tipotag(k).name ) ;
fprintf( fid,'%d %f %16.9e %f %f \n', xa(k), obj.f15.tipotag(k).val(1:2),obj.f15.tipotag(k).val(4:5) ) ;
end
fprintf( fid, '%d %s \n',  obj.f15.ntif , 'nbfr') ; %

for k = 1: obj.f15.ntif
lastr={obj.f15.tipotag.name};
xa(k)=str2double(lastr{1,k}(end));
xa(isnan(xa))=0; xa(xa>2)=0;
fprintf( fid, '%s \n',  obj.f15.tipotag(k).name ) ;
fprintf( fid,'%16.9e %f %f \n',obj.f15.tipotag(k).val(2),obj.f15.tipotag(k).val(4:5) ) ;
end

% fprintf(fid,'%d %s \n', size(obj.p,2)-1,'nope');
fprintf(fid,'%d %s \n', length(bounodes),'nope');

for hh=1:length(bounodes) % loop through all boundaries.
fprintf(fid,'%d %d %d %d %d %s %d %s  \n', length(bounodes{hh}), ops, '!', hh,'Open boundary');
% 
   
for k = 1: obj.f15.ntif
const{k}=lower(obj.f15.tipotag(k).name);
end 
pts_data=bounodes{hh}; % Number of boundary/nodes
[amp,phase ]= gen_harm_FES_SCHISM(iflag, const,pts_data);

  
for i=1:length(const)
  fprintf(fid,'%s\n',obj.f15.tipotag(i).name);
  fprintf(fid,'%f %f\n',[amp(:,i,1) phase(:,i,1)]');
end %for

if(iflag==2)
  for i=1:length(const)
    fprintf(fid,'%s\n',obj.f15.tipotag(i).name);
    fprintf(fid,'%f %f %f %f\n',[amp(:,i,2) phase(:,i,2) amp(:,i,3) phase(:,i,3)]');
  end %for
end
end 
if ops(3)==4
    fprintf(fid,'%s\n','1.0 ! T nudge');
end 

if ops(4)==4
    fprintf(fid,'%s\n','1.0 ! S nudge');
end 

fclose(fid); close all; 


function obj = tide_fac(ll,t_s,t_e,incnstit)
% obj = tide_fac(obj,t_s,t_e,incnstit)
% Input a msh class obj, start and time simulation times, and constituent
% names and return populated f15 tidal potential options in the msh class
% obj                
%                                                                       
% Created by William Pringle March 15 2018, 
% Copied functions from U_tide to calculate the tide factors, available at: 
% https://www.mathworks.com/matlabcentral/fileexchange/
% 46523--utide--unified-tidal-analysis-and-prediction-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the times and lat entries
% ts = datenum(t_s);
% te = datenum(t_e);
ts = datenum(t_s);
te = datenum(t_e);
tref = 0.5*(ts + te);
% make a fake t with 1000 intervals
t = linspace(ts,te,1000)';
lat = mean(ll(:,2));

%% Allow for an entry just to be major8 to get the major eight constituents
if strcmp(incnstit{1},'major8')
   incnstit = {'M2','S2','N2','K2','K1','O1','Q1','P1'};
   disp('Using the major eight harmonic tidal constituents')
   disp(incnstit)
end

%% Now check and get the indices of the incnstits (this part came from U_tide)
load('tide_fac_constants.mat','const');
cnstit.NR.lind = nan*ones(size(incnstit,1),1);
for j = 1:length(incnstit)
    lind1 = strmatch(incnstit{j},const.name);
    if isempty(lind1)
        error(['tide_fac: unrecognized non-reference constituent: '...
               incnstit{j} '.']);
    else
        cnstit.NR.lind(j) = lind1;
    end
end
[~,seq] = sort(const.freq(cnstit.NR.lind));
cnstit.NR.lind = cnstit.NR.lind(seq);
if isempty(cnstit.NR.lind)
    error('tide_fac: no constituents specified');
end
cnstit.NR.frq = const.freq(cnstit.NR.lind);
cnstit.NR.name = cellstr(const.name(cnstit.NR.lind,:));
% Doodson amp multiplied by Doodson constant/g:  
% (See : https://www.whoi.edu/fileserver.do?id=21351&pt=10&p=17272)
cnstit.NR.amp  = abs(const.doodsonamp(cnstit.NR.lind))*0.2675; 
cnstit.NR.red  = const.earthreduc(cnstit.NR.lind); 

%% Now get the F, U, and V (using the U_tide function ut_FUV)
[F,U,V] = ut_FUV(t,tref,cnstit.NR.lind,lat,zeros(1,4));

%% Final factors
F = mean(F)'; % This is the average nodal factor
phs = (mean(U)+V(1,:))'*360;  % This is the average nodal correction + 
                              % astronomical argument at beginning of
                              % simulation
% Make sure phase between 0 and 360 for aesthetic purposes
while any(phs < 0)
    phs(phs<0) = phs(phs<0) + 360;
end

obj.f15.ntif = length(cnstit.NR.name);
for k = 1:obj.f15.ntif
    obj.f15.tipotag(k).name   = cnstit.NR.name{k}; % name in frq order
    obj.f15.tipotag(k).val(1) = cnstit.NR.amp(k); % potential amplitude of species
    obj.f15.tipotag(k).val(2) = cnstit.NR.frq(k)*2*pi/3600; % frq in rad/s format
    obj.f15.tipotag(k).val(3) = cnstit.NR.red(k); % earth rigidity reduction factor
    obj.f15.tipotag(k).val(4) = F(k);   % average nodal factor
    obj.f15.tipotag(k).val(5) = phs(k); % average nodal correction + 
                         % astronomical argument at beginning of simulation
end
%EOF
end

% Stuff below has been copied and slightly adjusted from U_tide 
% (which was copied from T_tide)
function [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs)
% UT_FUV()
% compute nodal/satellite correction factors and astronomical argument
% inputs
%   t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%   lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
% output
%   F = real nodsat correction to amplitude [unitless] (nt x nc)
%   U = nodsat correction to phase [cycles] (nt x nc)
%   V = astronomical argument [cycles] (nt x nc)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)

nt = length(t);
nc = length(lind);
%% nodsat
if ngflgs(2) % none
    F = ones(nt,nc);
    U = zeros(nt,nc);
else
    if ngflgs(1) % linearized times
        tt = tref;
    else         % exact times
        tt = t;
    end
    ntt = length(tt);
    load('tide_fac_constants.mat');
    [astro,~]=ut_astron(tt');
    if abs(lat) < 0.5 % to make sure no close to zero values
        lat = sign(lat)*0.5; 
        if sign(lat) == 0
            lat = 0.5;
        end
    end
    slat=sin(pi*lat/180);
    rr=sat.amprat;
    j=find(sat.ilatfac==1);
    rr(j)=rr(j).*0.36309.*(1.0-5.0.*slat.*slat)./slat;
    j=find(sat.ilatfac==2);
    rr(j)=rr(j).*2.59808.*slat; 
    uu=rem( sat.deldood*astro(4:6,:)+sat.phcorr(:,ones(1,ntt)), 1);
    nfreq=length(const.isat); %#ok
    mat = rr(:,ones(1,ntt)).*exp(1i*2*pi*uu);
    F = ones(nfreq,ntt);
    ind = unique(sat.iconst);
    for i = 1:length(ind)
        F(ind(i),:) = 1+sum(mat(sat.iconst==ind(i),:),1);
    end
    U = imag(log(F))/(2*pi); % faster than angle(F)
    F=abs(F);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        exp2 = abs(exp1);
        F(k,:)=prod(F(j,:).^exp2(:,ones(ntt,1)),1);
        U(k,:)=sum(U(j,:).*exp1(:,ones(ntt,1)),1);
    end
    F=F(lind,:)';
    U=U(lind,:)';
    if ngflgs(1) % nodal/satellite with linearized times
        F = F(ones(nt,1),:);
        U = U(ones(nt,1),:);
    end
end
%% gwch (astron arg)
if ngflgs(4) % none (raw phase lags not greenwich phase lags)
    if ~exist('const','var')
        load('tide_fac_constants.mat','const');
    end
    [~,ader] = ut_astron(tref);
    ii=isfinite(const.ishallow); 
    const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
    for k=find(ii)'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
    end
    V = 24*(t-tref)*const.freq(lind)';
else 
    if ngflgs(3)  % linearized times
        tt = tref;
    else 
        tt = t;   % exact times
    end
    ntt = length(tt);
    if exist('astro','var')
        if ~isequal(size(astro,2),ntt)
            [astro,~]=ut_astron(tt');
        end        
    else
        [astro,~]=ut_astron(tt');
    end
    if ~exist('const','var')
        load('tide_fac_constants.mat');
    end
    V=rem( const.doodson*astro+const.semi(:,ones(1,ntt)), 1);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        V(k,:) = sum(V(j,:).*exp1(:,ones(ntt,1)),1);
    end
    V=V(lind,:)';
    if ngflgs(3)    % linearized times
        [~,ader] = ut_astron(tref);
        ii=isfinite(const.ishallow);
        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
        for k=find(ii)'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
                shallow.coef(ik) );
        end
        V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
    end
end
%EOF
end

function [astro,ader] = ut_astron(jd)
% UT_ASTRON()
% calculate astronomical constants
% input
%   jd = time [datenum UTC] (1 x nt)
% outputs
%   astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
%   ader = matrix of derivatives of astro [cycles/day] (6 x nt)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (copy of t_astron.m from t_tide, Pawlowicz et al 2002)

d=jd(:)'-datenum(1899,12,31,12,0,0);
D=d/10000;
args=[ones(size(jd));
      d;
      D.*D;
      D.^3];
sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];
astro=rem( [sc;hc;pc;npc;ppc]*args./360.0 ,1);
tau=rem(jd(:)',1)+astro(2,:)-astro(1,:);
astro=[tau;astro];
dargs=[zeros(size(jd));
       ones(size(jd));
       2.0e-4.*D;
       3.0e-4.*D.*D];
ader=[sc;hc;pc;npc;ppc]*dargs./360.0;
dtau=1.0+ader(2,:)-ader(1,:);
ader=[dtau;ader];
end


%% Read file hgrid.gr3 & set boundary flags
function [open_LL] = extract_open_nodes(filein)
% clear all
[EToV,VX,B,opedat,boudat,title] =readfort14(filein,1); % hgrid.ll  hgrid.gr3

ns=opedat.nvdll;
dta=opedat.nbdv;

for h=1:length(opedat.nbdv(1,:))
dta1{h}=dta(1:ns(h),h);
open_LL{h}=VX(dta1{h},:);
end 
clearvars -except open_LL
end 
%% READ GRID (*.GR3, FORT.14)
function [EToV,VX,B,opedat,boudat,title] = readfort14( finame, read_bou)
% Read ADCIRC grid file fort.14

disp('reading hgrid.gr3 file') ;
tic
if ( nargin == 0 )
    finputname = 'fort.14';
else
    finputname = finame ;
end

fid = fopen(finputname) ;

agrid = fgetl(fid) ;
disp(agrid) ;
title = agrid ;

N = fscanf(fid,'%g %g',2) ;

% Val = zeros(N(2),4) ;

%
% for i = 1: N(2)
%    Val(i,1:4) = fscanf(fid,'%d %g %g %g \n', 4 ) ;
% end
%
% Nov 15, 2012, improve reading efficiency
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;

%Val = Val(iv,:) ;

% idx = zeros(N(1),5) ;
%
% for i = 1: N(1)
%    idx(i,1:5) = fscanf(fid,'%d %d %d %d %d \n', 5) ;
% end
%
% Nov 15, 2012, improve reading efficient
idx = fscanf(fid,'%d %d %d %d %d \n', [5 N(1)])' ;

%idx = idx(iv,:) ;

% VX = zeros(N(2),2) ;
% B = zeros(N(2),2) ;
% EToV = zeros(N(1),3) ;

% Arrange it to a Nodal DG input
VX = Val(:,2:3) ;
B  = Val(:,4) ;
EToV = idx(:,3:5) ;

if(read_bou)
    % Read in boundary
    % Open boundary
    msgline = fgetl(fid) ;
    nope = sscanf(msgline,'%d %*s') ;
    
    msgline = fgetl(fid) ;
    neta = sscanf(msgline,'%d %*s') ;
    
    nvdll = zeros(nope,1) ;
    ibtypee = zeros(nope,1) ;
    nbdv = zeros(neta,nope) ;
    % tic
    for i = 1: nope
        msgline = fgetl(fid) ;
        
        [varg] = sscanf(msgline,'%d %*s \n') ;
        nvdll(i) = varg ;
        ibtypee(i) = 0 ;
        
        %
        % for k = 1: nvdll(i)
        %
        %    % % nbdv(k,i) = fscanf(fid,'%d \n')
        %    % % fscanf(fid,'%d \n')
        %    msgline = fgetl(fid) ;
        %
        %    nbdv(k,i) = str2num(msgline) ;
        % end
        %
        % Nov 25, 2012, improve reading efficiency
        nbdv(1:nvdll(i),i) = fscanf(fid,'%d \n', nvdll(i) ) ;
    end
    % toc
    
    % ocean boundary
    opedat.nope = nope ;
    opedat.neta = neta ;
    opedat.nvdll = nvdll ;
    opedat.ibtypee = ibtypee ;
    opedat.nbdv = nbdv ;
    
    % land boundary
    msgline = fgetl(fid) ;
    nbou = sscanf(msgline,'%d %*s') ;
    
    msgline = fgetl(fid) ;
    nvel = sscanf(msgline,'%d %*s') ;
    
    nvell = zeros(nbou,1) ;
    ibtype = zeros(nbou,1) ;
    nbvv = sparse(nvel,nbou) ;
    ibconn = sparse(nvel,nbou) ;
    barinht = sparse(nvel,nbou) ;
    barincfsb = sparse(nvel,nbou) ;
    barincfsp = sparse(nvel,nbou) ;
    
    % tic
    for i = 1: nbou
        msgline = fgetl(fid) ;
        
        [varg,~] = sscanf(msgline,'%d %d %*s \n') ;
        nvell(i) = varg(1) ;
        ibtype(i) = varg(2) ;
        
        switch ( ibtype(i) )
            case {0,1,2,10,11,12,20,21,22,30,60,61,101,52}
                % Nov 15, 2012, improve reading efficiency
                nbvv(1:nvell(i),i) = fscanf(fid,'%d \n', nvell(i) ) ;
            case  {3, 13, 23}
                disp('3 13 23')
                val = fscanf(fid,'%g %g %g \n', [3 nvell(i)] )  ;
                nbvv(1:nvell(i),i) = val(1,:) ;
            case  {4, 24}
                %disp('4 24')
                val = fscanf(fid,'%g %g %g %g %g \n', [5 nvell(i)] )  ;
                nbvv(1:nvell(i),i) = val(1,:) ;
                ibconn(1:nvell(i),i) = val(2,:) ;
                barinht(1:nvell(i),i) = val(3,:) ;
                barincfsb(1:nvell(i),i) = val(4,:) ;
                barincfsp(1:nvell(i),i) = val(5,:) ;
            case  {5, 25}
                %disp('5 25')
                val = fscanf(fid,'%g % g %g %g %g %g %g %g \n', [8 nvell(i)] ) ;
                nbvv(1:nvell(i),i) = val(1,:) ;
                %otherwise
                %    msgline = fgetl(fid) ;
        end
    end
    % toc
    
    % land boundary
    boudat.nbou = nbou ;
    boudat.nvel = nvel ;
    boudat.nvell = nvell ;
    boudat.ibtype = ibtype ;
    boudat.nbvv = nbvv ;
    
    if ( sum(ibtype == 24) > 0 ||  sum(ibtype == 4) > 0 )
        boudat.ibconn = ibconn ;
        boudat.barinht = barinht ;
        boudat.barincfsb = barincfsb ;
        boudat.barincfsp = barincfsp ;
    end
else
    opedat = []; 
    boudat  = []; 
end

fclose(fid) ;

return

end




function [amp,phase ]= gen_harm_FES_SCHISM(iflag, const,pts_data)
% iflag=2;
%Generate ampl. and phases for elev, u,v from FES2014 using linear interp directly
%from nc files.
%Requires inputs: (1) open.ll (generated by gen_fg.f90 with 1st 2 lines removed (to easily read into mlab), 
%                 1 seg at a time);  
%                 (2) *.nc (in this dir)
%                 (3) iflag: 1: generate elev only; 2: elev, u,v
%  Also need to first download FES2014 elev and velocity data at AVISO site
%Outputs: ap.out (for bctides.in; see below for order)

% open=load('fg.bp'); %ID,lon,lat of open bnd nodes
% npt=size(open,1);

npt=size(pts_data,1);

%miss_value=1.e10; %junk value
for ifl=1:2*iflag-1 %loop over elev, u,v
%-----------------------------------------------------------
for i=1:length(const)
  if(ifl==1)
    ncid = netcdf.open(['.\ocean_tide_extrapolated\' const{i} '.nc'],'NC_NOWRITE');
  elseif(ifl==2)
    ncid = netcdf.open(['.\eastward_velocity\' const{i} '.nc'],'NC_NOWRITE');
  elseif(ifl==3)
    ncid = netcdf.open(['.\northward_velocity\' const{i} '.nc'],'NC_NOWRITE');
  else
    error('Unknown iflag')
  end
  vid=netcdf.inqVarID(ncid,'lat'); %1D
  lat = netcdf.getVar(ncid, vid); %ascending order (-90,90)
  vid=netcdf.inqVarID(ncid,'lon');
  lon = netcdf.getVar(ncid, vid); %ascending order [0,360)
  nx=length(lon); ny=length(lat);

  if(ifl==1)
    vid=netcdf.inqVarID(ncid,'amplitude'); %(nx,ny)
  elseif(ifl==2)
    vid=netcdf.inqVarID(ncid,'Ua'); %(nx,ny)
  else
    vid=netcdf.inqVarID(ncid,'Va'); %(nx,ny)
  end
  amp0= netcdf.getVar(ncid, vid);
  amp0=amp0/100; %to meters or m/s

  if(ifl==1)
    vid=netcdf.inqVarID(ncid,'phase');
  elseif(ifl==2)
    vid=netcdf.inqVarID(ncid,'Ug');
  else
    vid=netcdf.inqVarID(ncid,'Vg');
  end
  pha0= netcdf.getVar(ncid, vid); %degr

  %Put phases into [0,360)
  indx=find(pha0<0);
  pha0(indx)=pha0(indx)+360;
  indx=find((pha0<0 | pha0>360) & pha0<1.4e10); %not junk but out of range
  if(sum(indx) ~=0) 
    pha0(indx)
    error('Phase out of bound');
  end

  ap=zeros(npt,2);
  for j=1:npt
%     lon2=open(j,2); if(lon2<0); lon2=lon2+360; end;
%     lat2=open(j,3);
    
    lon2=pts_data(j,1); if(lon2<0); lon2=lon2+360; end;
    lat2=pts_data(j,2);
    
    I=find(lon>=lon2,1); %find 1st entry
    J=find(lat>=lat2,1);
    if(I<=1 || J<=1) 
      disp('Failed to find an interval'); [j lon2 lat2 I J]
      error('Bomb out'); 
    end

    ratx=(lon(I)-lon2)/(lon(I)-lon(I-1));
    raty=(lat(J)-lat2)/(lat(J)-lat(J-1));
    %Check junks
    amp_max=max([amp0(I-1,J-1);amp0(I,J-1);amp0(I-1,J);amp0(I,J)]);
    amp_min=min([amp0(I-1,J-1);amp0(I,J-1);amp0(I-1,J);amp0(I,J)]);
    if(amp_max>100)
      if(amp_min>100 || amp_min<0)
        [i j I J]
        error('All junks for amp:'); 
      else %use min
        amp(j,i,ifl)=amp_min;
      end
    else
      tmp1=ratx*amp0(I-1,J-1)+(1-ratx)*amp0(I,J-1);
      tmp2=ratx*amp0(I-1,J)+(1-ratx)*amp0(I,J);
      amp(j,i,ifl)=tmp1*raty+tmp2*(1-raty);
    end %amp_max

%    if(ap(i,1)>2)
%      disp('Amp>2m');
%      [i j amp_max amp_min]
%    end

    %Make sure there is no jump in phase range
    wild(1:4)=[pha0(I-1,J-1); pha0(I,J-1); pha0(I,J); pha0(I-1,J)]; %counter-clockwise
    pha_max=max(wild);
    pha_min=min(wild);

    if(pha_max>370) %junk
      if(pha_min<0 ||  pha_min>370)
        [i j I J]
        error('All junk phases:'); 
      end
      phase(j,i,ifl)=pha_min;
    else
      for k=1:4
        if(abs(wild(k)-pha_max)>180)
          wild(k)=wild(k)+360;
          if(abs(wild(k)-pha_max)>180)
            [i j k I J]
            wild(1:4)
            error('Phases still jump');
          end
        end
      end %for k

      tmp1=ratx*wild(1)+(1-ratx)*wild(2);
      tmp2=ratx*wild(4)+(1-ratx)*wild(3);
      phase(j,i,ifl)=tmp1*raty+tmp2*(1-raty);
    end %pha_max

  end %for j=1:npt
end %for i - freq's
%-----------------------------------------------------------
end %for ifl

end
