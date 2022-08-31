%clear; close all;

%-------------------inputs---------------------------------------
%rundir with source_sink.in and vsource.th
wDir='../';

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

