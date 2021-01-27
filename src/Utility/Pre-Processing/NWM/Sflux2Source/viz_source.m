clear; close all;

%-------------------inputs---------------------------------------
%rundir with source_sink.in and vsource.th
wDir='/sciclone/data10/wangzg/NWM/RUN1i_ZG/';

hgrid_name='hgrid.ll'; %hgrid.ll or hgrid.gr3, must have bnd info
n_largest = 30; %mark the largest 'n_largest' sources and sinks; 
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

%read hgrid
[ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(wDir,hgrid_name,hgrid_name,2);
%calculate element center
node=[node; [nan nan nan nan]]; % add a dummy node at the end
ele(isnan(ele))=np+1;
x=nanmean(reshape(node(ele',2),4,ne)',2);
y=nanmean(reshape(node(ele',3),4,ne)',2);
node=node(1:end-1,:); ele(ele==np+1)=nan; %restore node and ele

ns=length(sid);

tt=vs(:,1)/86400;
% ms0=ms;
% ms=ms(:,2:ns+1);
vs=vs(:,2:end);
if ne_sink>0
    vi=vi(:,2:end);
end

%plot sources (and sinks if any) as patches
vs_ele=nan(ne,1); vs_ele(sid)=mean(vs,1);
patch('Faces',ele,'Vertices',node(:,2:3),'FaceVertexCData',vs_ele,'FaceColor','flat','LineStyle','none');
hold on;
colormap jet; caxis([0 1]); colorbar;

% scatter(x(sid(1:n_skip:end)), y(sid(1:n_skip:end)), max(0.001,mean(vs(:,1:n_skip:end),1)),'filled'); hold on;
% if (ne_sink > 0)
%     scatter(x(iid(1:n_skip:end)), y(iid(1:n_skip:end)), max(0.001,-mean(vi(:,1:n_skip:end),1))','filled'); hold on;
% end

%mark the largest n sources/sinks with a circle; the size of the circle is
%proportional to the time-average source/sink
[dummy, ilarge]=sort(mean(vs,1),'descend');
if (ne_sink > 0)
    [dummy, ivi]=sort(mean(abs(vi),1),'descend');
end
for i=1:n_largest
    scatter(x(sid(ilarge(i))),y(sid(ilarge(i))), max(0.001,mean(vs(:,ilarge(i)),1)),'MarkerEdgeColor',[1 0 0]); hold on;
    scatter(x(sid(ilarge(i))),y(sid(ilarge(i))),50,'+','MarkerEdgeColor',[1 0 0]); hold on;
    if (ne_sink > 0)
        scatter(x(iid(ivi(i))),y(iid(ivi(i))),50,'+','MarkerEdgeColor',[0 0 1]); hold on;
    end
end

axis equal;

