%clear; close all;

%-------------------inputs---------------------------------------
%rundir with source_sink.in and vsource.th
wDir=hdir;

hgrid_name='hgrid.ll'; %hgrid.ll or hgrid.gr3, must have bnd info
n_largest = 2; %mark the largest 'n_largest' sources and sinks; set to 0 to disable the marking
%-------------------end inputs---------------------------------------


% ms=load([wDir 'msource.th']);

%read source ele
[ne_source]=textread([wDir 'source_sink.in'], '%d', 1, 'headerlines', 0);
[sid]=textread([wDir 'source_sink.in'], '%d', ne_source, 'headerlines', 1);
vs=load([wDir 'vsource.th']);

%read sink ele if any
[ne_sink]=textread([wDir 'source_sink.in'], '%d', 1, 'headerlines', 1+ne_source+1);
if ne_sink>0
    [iid]=textread([wDir 'source_sink.in'], '%d', ne_sink, 'headerlines',1+ne_source+1+1);
    vi=load([wDir 'vsink.th']);
end

%read hgrid and save a copy
if exist([wDir 'hgrid.mat'],'file')
    load([wDir 'hgrid.mat']);%plot grid boundary
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
    [ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(wDir,hgrid_name,hgrid_name,2);
    save([wDir 'hgrid.mat'],'ne', 'np', 'node', 'ele', 'i34', 'bndnode','open_bnds','land_bnds','ilb_island');
end
x=nan(ne,1); y=x;
for i=1:ne
    x(i)=mean(node(ele(i,1:i34(i)),2));
    y(i)=mean(node(ele(i,1:i34(i)),3));
end

ns=length(sid);

tt=vs(:,1)/86400;
% ms0=ms;
% ms=ms(:,2:ns+1);
vs=vs(:,2:end);
if ne_sink>0
    vi=vi(:,2:end);
end

%plot sources (and sinks if any) as circles of different sizes
%a minimum size is forced.
scatter(x(sid), y(sid), max(0.001,mean(vs,1)),'filled'); hold on;
if (ne_sink > 0)
    scatter(x(iid), y(iid), max(0.001,-mean(vi,1))','filled'); hold on;
end

%mark the largest n sources/sinks with a fixed-size circle
[dummy, ivs]=sort(mean(vs,1),'descend');
if (ne_sink > 0)
    [dummy, ivi]=sort(mean(abs(vi),1),'descend');
end
for i=1:n_largest
    scatter(x(sid(ivs(i))),y(sid(ivs(i))),50,'+'); hold on;
    if (ne_sink > 0)
        scatter(x(iid(ivi(i))),y(iid(ivi(i))),50,'+'); hold on;
    end
end

axis equal;
















% figure;
% tt=[1:size(vs,1)]/24;
% for i=1:2
%     plot(tt,vs(:,ivs(i))); hold on;
% end




% II=find(sid==159);
% plot(tt,ms(:,II)); hold on;
% mean(vs(:,II))

% mean_vs=mean(vs,1);
% iRiver=find(mean_vs>0) %river inputs
% for k=1:20:365
%     clf;
%     scatter(x(sid(:)), y(sid(:)), max(vs(k,:),10), ms(k,:),'filled');
%     title(['Day ' num2str(tt(k)) '; ' datestr(datenum('2004-1-1')+tt(k))]);
%     caxis([0 30]);
%     colorbar;
% end
% figure
% for i=1:length(II)
%     plot(tt,ms(:,II(i))); hold on;
% end







% 
% %add sources from flux.th
% flux=load([wDir 'flux.th']);
% [n_flux]=textread([wDir 'flux.bp'], '%d', 1, 'headerlines', 1);
% [dummy, x_flux, y_flux , ~]=textread([wDir 'flux.bp'], '%d %f %f %f', n_flux, 'headerlines', 2);
% hold on; scatter(x_flux, y_flux, 10,'filled','r');
% 
% ie_flux=nan(n_flux,1);
% for i=1:n_flux
%     dist=((x-x_flux(i)).^2 + (y-y_flux(i)).^2).^0.5; 
%     [minDist ie_flux(i)]=min(dist);
%     II=find(sid==ie_flux(i));
%     if isempty(II)
%         sid=[sid; ie_flux(i)];
%     else
%         %add vsource and print out msource to check
%     end
% end
% 
% 
% 
% 
% 
% dlmwrite([wDir 'vsource.xyz'],[x(sid(iRiver)) y(sid(iRiver)) mean(vs(:,iRiver))'],'delimiter',' ','precision',15);
% 
% for i=1:length(iRiver)
%     dist(i,:)=( (xlb-xlb(iRiver(i))).^2+(ylb-ylb(iRiver(i))).^2 ).^0.5;
% end
% [dummy II]=min(dist,[],1);
% iSurro=iRiver(II);
% ms_surro=ms(:,iSurro);
% 
% % colormap('jet'); 
% % for k=1:20:365
% %     clf;
% % %     scatter(x(sid(:)), y(sid(:)), max(vs(k,:),3), ms(k,:),'filled'); hold on;
% %     scatter(x(sid(:)), y(sid(:)), max(vs(k,:),3), ms_surro(k,:),'filled'); hold on;
% %     scatter(x(sid(iRiver)),y(sid(iRiver)),'r+');
% %     title(['Day ' num2str(tt(k)) '; ' datestr(datenum('2004-1-1')+tt(k))]);
% %     caxis([0 30]);
% %     colorbar;
% % end
% 
% ms0(:,2:ns+1)=ms_surro;
% dlmwrite([wDir 'msource_surro.th'],ms0,'delimiter',' ','precision',15);
