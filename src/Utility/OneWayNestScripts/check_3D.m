%Plot uv3D.th as quiver, or time series  plots for T,S
%Inputs: fort.2[3-5] from interpolate_variables_selfe*.f90
clear all; close all;
xy=load('fort.25'); %ifile,nvrt_fg,xfg(),yfg()
u=load('fort.23'); %time (s), u|T@1,2,... (all levels)
v=load('fort.24'); %time (s), v|S@1,2,...

ifile=xy(1,1); %2: T,S; 3: uv
nvrt=xy(1,2);
ntime=size(u,1);
npt=(size(u,2)-1)/nvrt
x=xy(3:2+npt);
y=xy(3+npt:2+2*npt);
xmin=min(x); xmax=max(x);
ymin=min(y); ymax=max(y);
%Location of vel. scale
xlabel=0.8*xmin+0.2*xmax;
ylabel=-0.1*ymin+1.1*ymax;
scale_factor1=500; %adjust vel. scale for surface
scale_factor2=5000; %for bottom

if(ifile==2) %T,S
  for i=1:npt
    figure(1);
    subplot(7,7,i);
    plot(u(:,1)/86400,u(:,2+(i-1)*nvrt:1+i*nvrt),'k');
    title(['T; station ' num2str(i)]);

    figure(2);
    subplot(7,7,i);
    plot(v(:,1)/86400,v(:,2+(i-1)*nvrt:1+i*nvrt),'k');
    title(['S; station ' num2str(i)]);
  end %for
else %uv
  avi_out1 = avifile('uv3D_surf.avi','FPS',5);
  avi_out2 = avifile('uv3D_bot.avi','FPS',5);
  for it=1:ntime
    tt=u(it,1)/86400; %days
    %u2(1:nvrt,1:npt)
    u2=reshape(u(it,2:end),nvrt,npt);
    v2=reshape(v(it,2:end),nvrt,npt);

    figure(1); %surface
    quiver([x xlabel],[y ylabel],[u2(nvrt,:) 0.5]*scale_factor1,[v2(nvrt,:) 0]*scale_factor1,'AutoScale','off');
    axis([1.2*xmin-0.2*xmax -0.2*xmin+1.2*xmax 1.2*ymin-0.2*ymax -0.2*ymin+1.2*ymax]);
    title(['t= ' num2str(tt) ' days; surface']);
    text(xlabel,ylabel,'0.5m/s');
    frame = getframe(gcf);
    avi_out1=addframe(avi_out1,frame);
    
    figure(2); %bottom
    quiver([x xlabel],[y ylabel],[u2(1,:) 0.1]*scale_factor2,[v2(1,:) 0]*scale_factor2,'AutoScale','off');
    axis([1.2*xmin-0.2*xmax -0.2*xmin+1.2*xmax 1.2*ymin-0.2*ymax -0.2*ymin+1.2*ymax]);
    title(['t= ' num2str(tt) ' days; bottom']);
    text(xlabel,ylabel,'0.1m/s');
    frame = getframe(gcf);
    avi_out2=addframe(avi_out2,frame);
    
    %return;
    figure(1); clf; 
    figure(2); clf;
   
  end %for it
  avi_out1=close(avi_out1);
  avi_out2=close(avi_out2);
end %ifile
