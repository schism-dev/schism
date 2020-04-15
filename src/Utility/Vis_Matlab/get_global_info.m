function [ne,np,nvrt,nm,xy00,i23d,ivs,vzcor,vid5,vdry_e,h0,dp,kbp00,vid_eta,nproc,np_lcl,ne_lcl,ns_lcl,iplg,ielg,iegl_rank]=get_global_info(icomb,base,varname,stacks)
%Get basic info before reading (combined or uncombined) nc outputs 
%Outputs starting from npoc are only defined if icomb=0 (uncombined)
%vzcor,vis5,vdry_e: IDs for zcor and varname, and wetdry_elem

if(icomb==0)
  fid=fopen([base '/outputs/local_to_global_0000'],'r');
  first=fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',18);
  %disp(first);
  fclose(fid);

  ns=first(1); ne=first(2); np=first(3); nvrt=first(4); nproc=first(5);
  nm=NaN(4,ne); i34=zeros(ne,1); xy00=zeros(np,2); iegl_rank=NaN(ne,1);

  for irank=0:nproc-1
    fid=fopen([base '/outputs/local_to_global_' num2str(irank,'%04.f')],'r');
    first=fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',18);
    %char=fscanf(fid,'%s\n')
    ne_lcl(irank+1)=fscanf(fid,'local to global mapping: %d\n',1);
    clear tmp2;
    tmp2=fscanf(fid,'%d %d',[2 ne_lcl(irank+1)]);
    ielg(irank+1,1:ne_lcl(irank+1))=tmp2(2,:);
    iegl_rank(tmp2(2,:))=irank;
    
    np_lcl(irank+1)=fscanf(fid,'%d\n',1);
    clear tmp2;
    tmp2=fscanf(fid,'%d %d',[2 np_lcl(irank+1)]);
    iplg(irank+1,1:np_lcl(irank+1))=tmp2(2,:);

    ns_lcl(irank+1)=fscanf(fid,'%d\n',1);
    clear tmp2;
    tmp2=fscanf(fid,'%d %d',[2 ns_lcl(irank+1)]);
    islg(irank+1,1:ns_lcl(irank+1))=tmp2(2,:);

    tmp=fscanf(fid,'%*s %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n',16);
    %tmp=fscanf(fid,'Header:')
    %tmp=fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n',16)
%    tmp=fscanf(fid,'%d %f %d %d %d %f %f %f %f %f %f\n',1)
    nrec=fix(tmp(6)); dtout=tmp(7); kz=fix(tmp(10)); h0=tmp(11);
    format='';
    for j=1:nvrt
      format=[format '%f'];
    end
    format=[format '\n'];
    tmp=fscanf(fid,format,nvrt);

    tmp=fscanf(fid,'%d %d\n',2);
    clear tmp2;
    tmp2=fscanf(fid,'%g %g %g %g',[4 np_lcl(irank+1)]);
    xy00(iplg(irank+1,1:np_lcl(irank+1)),1)=tmp2(1,:);
    xy00(iplg(irank+1,1:np_lcl(irank+1)),2)=tmp2(2,:);
    dp(iplg(irank+1,1:np_lcl(irank+1)))=tmp2(3,:);
    kbp00(iplg(irank+1,1:np_lcl(irank+1)))=fix(tmp2(4,:));

    clear tmp;
    for j=1:ne_lcl(irank+1)
      iegb=ielg(irank+1,j);
      i34(iegb)=fscanf(fid,'%d',1);
      for k=1:i34(iegb)
        tmp=fscanf(fid,'%d',1);
        nm(k,iegb)=iplg(irank+1,tmp);
      end %for k
    end %for j

    fclose(fid);
  end %for irank

  if(sum(isnan(iegl_rank))~=0)
    error(['get_global_info: iegl_rank has nan']);
  end

  %Read header of nc
  fname=[base,'/outputs/' 'schout_0000_' num2str(stacks(1)) '.nc'];
  ncid0 = netcdf.open(fname,'NC_NOWRITE');
  vid5=netcdf.inqVarID(ncid0,varname);
  i23d=netcdf.getAtt(ncid0,vid5,'i23d'); %1: 2D; 2: 3D (whole level)
  ivs=netcdf.getAtt(ncid0,vid5,'ivs');
  vzcor=netcdf.inqVarID(ncid0,'zcor');
  vdry_e=netcdf.inqVarID(ncid0,'wetdry_elem');
  vid_eta=netcdf.inqVarID(ncid0,'elev'); 
  netcdf.close(ncid0);

else %combined nc
  %Return junks for vars for icomb=0
  nproc=0; np_lcl=0; ne_lcl=0; ns_lcl=0; iplg=0; ielg=0; iegl_rank=0;

  fname=[base,'/outputs/' 'schout_' num2str(stacks(1)) '.nc'];
  ncid0 = netcdf.open(fname,'NC_NOWRITE');
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
  vid=netcdf.inqVarID(ncid0,'node_bottom_index');
  kbp00=double(netcdf.getVar(ncid0, vid)); %(np)

  vid5=netcdf.inqVarID(ncid0,varname); 
  i23d=netcdf.getAtt(ncid0,vid5,'i23d'); %1: 2D; 2: 3D (whole level)
  ivs=netcdf.getAtt(ncid0,vid5,'ivs');
  vzcor=netcdf.inqVarID(ncid0,'zcor'); 
  vdry_e=netcdf.inqVarID(ncid0,'wetdry_elem');
  vid_eta=netcdf.inqVarID(ncid0,'elev'); 
  netcdf.close(ncid0);

  %Pad conn table with NaN for tri's
  nm(find(nm<0))=nan;
  ne=size(nm,2);
end %icomb

end %function
