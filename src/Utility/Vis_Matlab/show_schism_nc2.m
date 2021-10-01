%Read SCHISM node-based outputs from scribe I/O versions and display
%Inputs: out2d*.nc and the corresponding nc output for the var
clear all; close all;
start_stack=6; end_stack=6;
ivs=1; %1: scalar; 2: vector
filename='out2d'; %nc file name
varname='sigWaveHeight'; %var name
levelout=1; %use 1 for 2D var
ispher_nowrap=0; %1: remove wrap around elem on the globe (jump can be any lon)
if(ivs==2) 
  filename2='';
  varname2='';
end

%xyz=load('hgrid.xyz'); %xyz part of .gr3
%ax=[2.e5  4.e5 1.e5 4.e5]; %axis appearance
%xcenter=(ax(1)+ax(2))/2;
%ycenter=(ax(3)+ax(4))/2;

scrsz = get(0,'ScreenSize');
%4 parameters of position: left bottom_coord width height
figure('Position',[1 scrsz(4)*0.2 scrsz(3)/2 scrsz(4)*0.7]);
% Read the variable and the vertical grid
delete('slab.avi');
vidObj = VideoWriter('slab.avi');
vidObj.FrameRate = 30;  % Default 30; smaller ->slower
%vidObj.Quality = 50;    % Default 75
open(vidObj);

for istack=start_stack:end_stack
  if(istack==start_stack)
%----------------------------------------------------
  ncid0 = netcdf.open(['outputs/out2d_' num2str(istack) '.nc'],'NC_NOWRITE');
  dimid = netcdf.inqDimID(ncid0,'nSCHISM_hgrid_node');
  [~,np] = netcdf.inqDim(ncid0,dimid);
  dimid = netcdf.inqDimID(ncid0,'nSCHISM_vgrid_layers');
  [~,nvrt] = netcdf.inqDim(ncid0,dimid);
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_face_nodes'); 
  nm=double(netcdf.getVar(ncid0, vid)); %(4,ne)
  ne=size(nm,2);
  vid=netcdf.inqVarID(ncid0,'time'); 
  time=double(netcdf.getVar(ncid0, vid));
  ntime=length(time);

  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_node_x'); 
  x=double(netcdf.getVar(ncid0, vid));
  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_node_y'); 
  y=double(netcdf.getVar(ncid0, vid));
%  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_edge_x'); 
%  xcj=double(netcdf.getVar(ncid0, vid));
%  vid=netcdf.inqVarID(ncid0,'SCHISM_hgrid_edge_y'); 
%  ycj=double(netcdf.getVar(ncid0, vid));
  netcdf.close(ncid0);

  %Deal with quads
  i34(1:ne)=nan;
  for i=1:ne
    if(nm(4,i)<0)
      i34(i)=3;
      nm(4,i)=nan;
    else
      i34(i)=4;
    end
  end %for

  %Take care of wrap around elem
  icolor_nd(1:np)=0;
  if(ispher_nowrap==1)
    %If a node is on dateline, make it 180 deg
    xtmp=x;
    xtmp(find(xtmp==-180))=180;
    x=xtmp;
    for i=1:ne
      lon_min=1.e10; lon_max=-1.e10;
      lon_min=min(x(nm(1:i34(i),i)));
      lon_max=max(x(nm(1:i34(i),i)));
      if(abs(lon_min-lon_max)>100)
        icolor_nd(nm(1:i34(i),i))=1;
      end
    end %for i
    uv=nan(np,2);
  end %ispher_nowrap==1
%----------------------------------------------------
  end %if 1st stack

  ncid4 = netcdf.open(['outputs/' filename '_' num2str(istack) '.nc'],'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid4,varname); %(*,np,ntime)
  i23d= netcdf.getAtt(ncid4,vid,'i23d');
  if(ivs==2)
    ncid5 = netcdf.open(['outputs/' filename2 '_' num2str(istack) '.nc'],'NC_NOWRITE');
    vid5=netcdf.inqVarID(ncid5,varname2); %(*,np,ntime)
  end

  for it=1 %ntime
    if(i23d==1) %2D
      uv(:,1)=double(netcdf.getVar(ncid4, vid,[0 it-1],[np 1]));
      if(ivs==2)
        uv(:,2)=double(netcdf.getVar(ncid5, vid5,[0 it-1],[np 1]));
      end
    else %3D
      uv(:,1)=double(netcdf.getVar(ncid4, vid,[levelout-1 0 it-1],[1 np 1]));
      if(ivs==2)
        uv(:,2)=double(netcdf.getVar(ncid5, vid5,[levelout-1 0 it-1],[1 np 1]));
      end
    end

    %Deal with junks
    uv(find(abs(uv)>1.e10))=nan;

    if(ivs==2) %vector
      scale=2e2; %scale to fit
      quiver(x,y,squeeze(uv(:,1))*scale,squeeze(uv(:,2))*scale,0,'b');
%      text(3.5e5,4.e5,'1 m/s');
    else %scalar
      surf=squeeze(uv(:,1)); 
      surf(find(icolor_nd==1))=nan;

      h1=patch('Faces',nm(:,:)','Vertices',[x y],'FaceVertexCData',surf,'FaceColor','interp','EdgeColor','none');
      caxis([0 10]);
      colormap(jet(40));
      colorbar;
    end %%scalar|vetor
%   axis(ax);
    title([varname '; Time=' num2str(time(it))]);

    % Add image to avi file
    set(gcf,'Color',[1 1 1]);
    set(gcf,'nextplot','replacechildren');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
%    clf; %clear figure to avoid overlay

  end %for it
  netcdf.close(ncid4);
end %for istack
close(vidObj);
