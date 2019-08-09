%Animation time series of elev. to go together with ArcMap animations
clear all; close all;
eta=load('RUN09b_elev_anim.dat'); %time (days), eta@1,2,....
scrsz = get(0,'ScreenSize'); %screen size
figure('Position',[1 scrsz(4)/5 scrsz(3)/5 scrsz(4)/5]);
avif={'Siuslaw_mouth','Siuslaw_Hwy101','Ump_mouth','Ump_Hwy101'};

nsta=size(eta,2)-1;
nstep=size(eta,1);
if(nsta ~= length(avif)); error('# ofavi files do not match'); end;
for i=1:nsta
  avi_out = avifile([avif{i} '-ts.avi'],'fps',19); %20 fps is default in ArcMap
  
  clf;
  for j=1:nstep
    h=plot(eta(:,1)*24,eta(:,i+1),'k',eta(j,1)*24*ones(1,2),[-1.e4 1.e4],'r');
    axis([0 2 -10 20]);
%    xlabel('Time (hours)'); title('Elev (m MHHW)');
    set(h,'LineWidth',4);
    set(gca,'FontSize',25);
    set(gcf,'Color',[1 1 1]);    

    frame = getframe(gcf);
    avi_out=addframe(avi_out,frame)
  end %for j
  avi_out=close(avi_out);
end %for i
