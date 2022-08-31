function [file_nm,i23d,ne,np,nvrt,nm,xy00,vdry_e,h0,dp,kbp00,vid_eta]=get_global_info2(base,varname,stacks)
%Get basic info before reading scribed outputs.
%vzcor,vis5,vdry_e: IDs for zcor and varname, and wetdry_elem
%file_nm: file name that contains varname (e.g. 'out2d')
%xy00(np,2); nm(4,ne)

  fname0=[base,'/outputs/' 'out2d_' num2str(stacks(1)) '.nc'];
  ncid0 = netcdf.open(fname0,'NC_NOWRITE');
  dimid = netcdf.inqDimID(ncid0,'nSCHISM_hgrid_node');
  [~,np] = netcdf.inqDim(ncid0,dimid);
  dimid = netcdf.inqDimID(ncid0,'nSCHISM_vgrid_layers');
  [~,nvrt] = netcdf.inqDim(ncid0,dimid);
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_face_nodes');
  nm=double(netcdf.getVar(ncid0, vid)); %(4,ne)
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_node_x');
  xy00(1:np,1)=double(netcdf.getVar(ncid0, vid)); 
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_node_y');
  xy00(1:np,2)=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'minimum_depth');
  h0=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'depth');
  dp=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'bottom_index_node');
  kbp00=double(netcdf.getVar(ncid0, vid)); %(np)

  vdry_e=netcdf.inqVarID(ncid0,'dryFlagElement');
  vid_eta=netcdf.inqVarID(ncid0,'elevation'); 
  netcdf.close(ncid0);

%  fname=[base,'/outputs/' 'zCoordinates_' num2str(stacks(1)) '.nc'];
%  ncid0 = netcdf.open(fname,'NC_NOWRITE');
%  vzcor=netcdf.inqVarID(ncid0,'zCoordinates'); 
%  netcdf.close(ncid0);
%
  file_nm=varname;
  fname=[base,'/outputs/' file_nm '_' num2str(stacks(1)) '.nc'];
  if(~exist(fname)) %try out2d
    file_nm='out2d';
    fname=[base,'/outputs/out2d_' num2str(stacks(1)) '.nc'];
  end
  ncid_var = netcdf.open(fname,'NC_NOWRITE');
  vid5=netcdf.inqVarID(ncid_var,varname); 
  i23d=netcdf.getAtt(ncid_var,vid5,'i23d'); %1: 2D; 2: 3D (whole level)
%  ivs=netcdf.getAtt(ncid_var,vid5,'ivs');
  netcdf.close(ncid_var);

  %Pad conn table with NaN for tri's
  nm(find(nm<0))=nan;
  ne=size(nm,2);

end %function
