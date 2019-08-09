%Read SCHISM nc outputs to compare with vis6
clear all; close all;

ax=[324843.715  366235.556 277153.749 306866.564]; %axis appearance

xcenter=(ax(1)+ax(2))/2;
ycenter=(ax(3)+ax(4))/2;

%Test pre-comb
if(1==2)
  %local_id=355; %local node ID in a rank
  %ncid0 = netcdf.open('schout_0011_2.nc','NC_NOWRITE');
  local_id=166; %local node ID in a rank
  ncid0 = netcdf.open('schout_0022_2.nc','NC_NOWRITE');

  vid=netcdf.inqVarID(ncid0,'uveln'); %input var./array name
  u= double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'vveln'); %input var./array name
  v= double(netcdf.getVar(ncid0, vid));
  netcdf.close(ncid0);
  u(local_id,end,end)
  v(local_id,end,end)
end %1==2

%Test post-comb
if(1==1)
  ncid0 = netcdf.open('outputs/schout_3.nc','NC_NOWRITE');
  dimid = netcdf.inqDimID(ncid0,'nSCHISM_hgrid_node');
  [~,np] = netcdf.inqDim(ncid0,dimid);
  dimid = netcdf.inqDimID(ncid0,'nSCHISM_vgrid_layers');
  [~,nvrt] = netcdf.inqDim(ncid0,dimid);
  vid=netcdf.inqVarID(ncid0,'time'); 
  time=double(netcdf.getVar(ncid0, vid));
  ntime=length(time);
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_face_nodes'); 
  nm=double(netcdf.getVar(ncid0, vid)); %(4,ne)
  vid=netcdf.inqVarID(ncid0,'hvel'); %(2,nvrt,np,ntime)
  %Warning: start in getVaR follows C convention and starts from 0!
  uv=double(netcdf.getVar(ncid0, vid,[0 nvrt-1 0 ntime-1],[2 1 np 1]));
  vid=netcdf.inqVarID(ncid0,'salt'); %(nvrt,np,ntime)
  S=double(netcdf.getVar(ncid0, vid,[nvrt-1 0 ntime-1],[1 np 1]));
  vid=netcdf.inqVarID(ncid0,'temp'); 
  T=double(netcdf.getVar(ncid0, vid,[nvrt-1 0 ntime-1],[1 np 1]));
  %Deal with junks
  uv(find(abs(uv)>1.e10))=nan;
  T(find(abs(uv)>1.e10))=nan;
  S(find(abs(uv)>1.e10))=nan;

  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_node_x'); 
  x=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_node_y'); 
  y=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_edge_x'); 
  xcj=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_edge_y'); 
  ycj=double(netcdf.getVar(ncid0, vid));
%  vid=netcdf.inqVarID(ncid0,'dahv'); 
%  dahv=double(netcdf.getVar(ncid0, vid));
%  vid=netcdf.inqVarID(ncid0,'hvel_side'); 
%  suv=double(netcdf.getVar(ncid0, vid)); %(2,nvrt,ns,ntime)
%  size(suv) %query dims
%  vid=netcdf.inqVarID(ncid0,'salt_elem'); 
%  S_elem=double(netcdf.getVar(ncid0, vid));
%  vid=netcdf.inqVarID(ncid0,'ICE_tracer_1'); 
%  hice=double(netcdf.getVar(ncid0, vid));
  netcdf.close(ncid0);

%  figure(1); hold on;
%  patch('Faces',nm(1:3,:)','Vertices',[x y],'FaceVertexCData',hice(:,end),'FaceColor','interp','EdgeColor','none');
%  caxis([0 3]);
%  axis(ax);
%  colormap(jet(40));
%  colorbar;
%  title(['Ice volume; Time=' num2str(time(end)/86400)]);
%
%  return;

  %axis([1.e5 5e5 0.5e5 5e5]);
  figure(1); hold on;
  scale=2e3; %scale to fit
  quiver(x,y,squeeze(uv(1,1,:))*scale,squeeze(uv(2,1,:))*scale,0,'b');
%  quiver(xcj,ycj,squeeze(suv(1,end,:,end))*scale,squeeze(suv(2,end,:,end))*scale,0,'b');
  quiver(3.5e5,4.e5,1*scale,0.,0,'r');
  text(3.5e5,4.e5,'1 m/s');
  axis(ax);
  title(['suv; Time=' num2str(time(end))]);

%  figure(2); hold on;
%  surf=squeeze(S_elem(end,:,end));
%  patch('Faces',nm(1:3,:)','Vertices',[x y],'FaceVertexCData',surf','FaceColor','flat','EdgeColor','none');
%  caxis([0 30]);
%  axis(ax);
%  colormap(jet(40));
%  colorbar;
%  title(['S@elem; Time=' num2str(time(end))]);

  figure(3); hold on;
  surf=S; 
  patch('Faces',nm(1:3,:)','Vertices',[x y],'FaceVertexCData',surf','FaceColor','interp','EdgeColor','none');
  caxis([0 30]);
  axis(ax);
  colormap(jet(40));
  colorbar;
  title(['S; Time=' num2str(time(end))]);
end %1==2

%Pre-comb hotstart
if(2==1)
  ncid0 = netcdf.open('outputs/hotstart_0020_200.nc','NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); 
  time=double(netcdf.getVar(ncid0, vid))
  vid=netcdf.inqVarID(ncid0,'xcj'); 
  xcj=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'ycj'); 
  ycj=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'su2'); %Note the index order (nvrt,ns)
  su2=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'sv2'); 
  sv2=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'uu2'); %(nvrt,np)
  uu2=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'vv2'); 
  vv2=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'xnd'); 
  xnd=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'ynd'); 
  ynd=double(netcdf.getVar(ncid0, vid));

  netcdf.close(ncid0);

  scale=1e0; %scale to fit
%  quiver(xcj,ycj,su2(:,end)*scale,sv2(:,end)*scale); %,0,'b');
%  ax2=axis;
%  xcenter=(ax(1)+ax(2))/2;
%  ycenter=(ax(3)+ax(4))/2;
%  quiver(xcenter,ycenter,1*scale,0.,0,'r');
%  text(xcenter,ycenter+50,'1 m/s');
%  axis(ax2);
%  axis(ax);
%  title(['suv; Time=' num2str(time)]);

  fid=fopen('tmp.xyuv','w');
  fprintf(fid,'%f %f %f %f\n',[xcj ycj su2(end,:)' sv2(end,:)']');
  fclose(fid);
end %1==2

%Post combined hotstart
if(2==1)
  ncid0 = netcdf.open('outputs/hotstart_it=200.nc','NC_NOWRITE');
  %vid=netcdf.inqVarID(ncid0,'tr_el');
  %tr_el=double(netcdf.getVar(ncid0, vid)); %(ntracers,nvrt,ne)
  vid=netcdf.inqVarID(ncid0,'su2');
  su2=double(netcdf.getVar(ncid0, vid)); %(nvrt,ns)
  vid=netcdf.inqVarID(ncid0,'sv2');
  sv2=double(netcdf.getVar(ncid0, vid)); %(nvrt,ns)
  vid=netcdf.inqVarID(ncid0,'tr_nd');
  tr_nd=double(netcdf.getVar(ncid0, vid)); %(ntracers,nvrt,ns)
  vid=netcdf.inqVarID(ncid0,'dfv');
  dfv=double(netcdf.getVar(ncid0, vid)); %(nvrt,ns)
  netcdf.close(ncid0);

  side=load('side.xy');
  fid=fopen('tmp.xyuv','w');
  fprintf(fid,'%f %f %f %f\n',[side(:,2:3) su2(end,:)' sv2(end,:)']');
  fclose(fid);
%  return;

  fid=fopen('hgrid.gr3','r');
  char=fgetl(fid);
  tmp1=str2num(fgetl(fid));
  ne=fix(tmp1(1));
  np=fix(tmp1(2));
  bathy(1:np,1:4)=nan;
  nm(1:ne,1:5)=nan;
  for i=1:np
    tmp=str2num(fgetl(fid));
    bathy(i,1:3)=tmp(1:3);
  end %for i
  for i=1:ne
    nm(i,:)=str2num(fgetl(fid));
  end %for i
  fclose(fid);
  bathy(:,4)=squeeze(tr_nd(2,end,:));
  
  figure(4); hold on;
  %plot(bnd(:,1),bnd(:,2),'k.','MarkerSize',2);
  patch('Faces',nm(:,3:5),'Vertices',bathy(:,2:3),'FaceVertexCData',bathy(:,4),'FaceColor','interp','EdgeColor','none');
  %surf=squeeze(tr_el(2,end,:));
  %patch('Faces',nm(:,3:5),'Vertices',bathy(:,2:3),'FaceVertexCData',surf,'FaceColor','flat','EdgeColor','none');
  colormap(jet(40));
  caxis([0 30]);
  colorbar;
  axis(ax);
  set(gcf,'Color',[1 1 1]);

  %disp('sample elev');
  %eta2(10697:10699)

end %1==2
