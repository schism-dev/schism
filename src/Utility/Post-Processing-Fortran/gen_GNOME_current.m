%Reformat SCHISM .nc outputs to GNOME format
%Inputs: outputs/schout_*.nc; nbe.out & GNOME_bnd.out (from gen_GNOME_bnd.f90); hgrid.ll.mlab (just node part)
%        Also change # of stacks below, and may adjust starting time also.
%Outputs: schism_gnome_*.nc (only 1st 3 nodes of quads are used but quads are not tested)

close all; clear all;
stack_start=1;
stack_end=2;
nbe=load('nbe.out'); %(ne,4) - 1st column is elem #
bnd=transpose(load('GNOME_bnd.out')); %bnd(4,:)
nbnd=size(bnd,2);
ll=load('hgrid.ll.mlab'); %ID, lon,lat,depth

for ifile=stack_start:stack_end %stack #s for SCHISM outputs
%--------------------------------------------------
filen=['outputs/schout_' num2str(ifile) '.nc'];
disp(['doing ' filen]);
ncid0 = netcdf.open(filen,'NC_NOWRITE');
[dimname, np] = netcdf.inqDim(ncid0,0);
[dimname, ne] = netcdf.inqDim(ncid0,1);
vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
time= double(netcdf.getVar(ncid0, vid));
ntime=length(time);

vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_face_nodes');
elnode=netcdf.getVar(ncid0, vid); %(4,ne)
vid=netcdf.inqVarID(ncid0,'hvel'); %(2,nvrt,np,ntime)
uv= double(netcdf.getVar(ncid0, vid));
netcdf.close(ncid0);

%Output
ncid2 = netcdf.create(['schism_gnome_' num2str(ifile) '.nc'],'CLOBBER');
varid = netcdf.getConstant('GLOBAL');
%netcdf.putAtt(ncid2,varid,'file_type','FEM');
%netcdf.putAtt(ncid2,varid,'Conventions','COARDS');
netcdf.putAtt(ncid2,varid,'grid_type','Triangular');

dims(1)=netcdf.defDim(ncid2,'node',np);
dims(2)=netcdf.defDim(ncid2,'nele',ne);
dims(3)=netcdf.defDim(ncid2,'nbi',4);
dims(4)=netcdf.defDim(ncid2,'time',ntime); %netcdf.getConstant('NC_UNLIMITED'));
dims(5)=netcdf.defDim(ncid2,'nbnd',nbnd);
dims(6)=netcdf.defDim(ncid2,'three',3);

% Define new variable in the file.
timeid = netcdf.defVar(ncid2,'time','float',dims(4));
netcdf.putAtt(ncid2,timeid,'long_name','Time');
%Adjust starting time
netcdf.putAtt(ncid2,timeid,'units','days since 2002-04-30 0:00:00 00:00');
netcdf.putAtt(ncid2,timeid,'base_date',int32([2002 4 30 0]));
latid = netcdf.defVar(ncid2,'lat','float',dims(1));
lonid = netcdf.defVar(ncid2,'lon','float',dims(1));
uid = netcdf.defVar(ncid2,'u','float',dims([1 4]));
vid = netcdf.defVar(ncid2,'v','float',dims([1 4]));
bnd_id=netcdf.defVar(ncid2,'bnd','int',dims([3 5]));
table_id=netcdf.defVar(ncid2,'nv','int',dims([2 6]));
%netcdf.putAtt(ncid2,table_id,'order','ccw');
nbe_id=netcdf.defVar(ncid2,'nbe','int',dims([2 6]));
netcdf.putAtt(ncid2,nbe_id,'order','ccw');

% Leave define mode and enter data mode to write data.
netcdf.endDef(ncid2)

netcdf.putVar(ncid2,bnd_id,bnd);
netcdf.putVar(ncid2,table_id,elnode(1:3,:)');
netcdf.putVar(ncid2,nbe_id,nbe(:,2:4));
netcdf.putVar(ncid2,timeid,single(time/86400));
netcdf.putVar(ncid2,latid,single(ll(:,3)));
netcdf.putVar(ncid2,lonid,single(ll(:,2)));
netcdf.putVar(ncid2,uid,single(uv(1,end,:,:)));
netcdf.putVar(ncid2,vid,single(uv(2,end,:,:)));
netcdf.close(ncid2);
%--------------------------------------------------
end %for ifile
