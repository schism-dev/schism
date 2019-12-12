%--------------------------------------------------------------------------
% Script for directing freshwater from NWM segments to SCHISM
%   by Fei Ye (feiye@vims.edu), 2019
%--------------------------------------------------------------------------
%----- (1) inputs:-----
%hgrid dir
hdir='./';
%nwm dataset dir
REPODir='/sciclone/data10/feiye/DELAWARE_REPO/00_TWL_Shared/01_data/01-NWM-4-isabel-irene-sandy-13sep2018/NetCDF_HDF/Irene_output/';

% Time origin in GMT 
TimeStart=datenum('2011-7-27');
TimeEnd=datenum('2011-7-27')+80;

%SCHISM boundary segments
bnds={
[45995 281534];    %start/end of a portion of outer boundary
[281535 273573]    %start/end of a portion of island boundary
};
%-----end inputs-----


%susq. r. element'
ie_plus_ll={[-76.13 39.62]};
susq_r_flow = [1500]; %a constant flow imposed at Susquehanna River

skip_segs=[]; 


%-----(2) Read hgrid and build grid topography-----
if exist([hdir 'hgrid.mat'],'file')
    load([hdir 'hgrid.mat']);
    %plot grid boundary
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
    [ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(hdir,'hgrid.ll','hgrid.ll',2);
    save([hdir 'hgrid.mat'],'ne', 'np', 'node', 'ele', 'i34','bndnode','open_bnds','land_bnds','ilb_island');
end

%find additional point sources where constant streamflows are imposed
n_plus = length(ie_plus_ll); 
for k=1:n_plus
  ele_x=nan(ne,1); ele_y=nan(ne,1);
  for i=1:ne
      ele_x(i,1)=mean( node(ele(i,1:i34(i)) , 2 ) );
      ele_y(i,1)=mean( node(ele(i,1:i34(i)) , 3 ) );
  end
  [dummy ie_plus(k)]=min((ele_x-ie_plus_ll{k}(1)).^2 + (ele_y-ie_plus_ll{k}(2)).^2 ); 
end
ie_plus

%diagnostic plot
%find land bnd nodes that needs point sources 
bnd_node=[];
for k=1:length(bnds)
    for i=1:length(land_bnds)
        [II idx]=ismember(bnds{k},land_bnds{i});
        if (sum(II) == 2 )
            if ilb_island(i)==0 %not an island
                bnds_nd=land_bnds{i}(idx(1):idx(2));
            else
                if idx(1)>idx(2) %not passing the end of the circle
                    bnds_nd=flip(land_bnds{i}( [idx(2):-1:1 length(land_bnds{i}):-1:idx(1) ] ) );
                elseif idx(1)<idx(2)
                    bnds_nd=land_bnds{i}([idx(1):idx(2)]);
                end
            end
            bnd_node=[bnd_node; bnds_nd];
            break;
        end %if found
    end %i
end %k
line(node(bnd_node,2),node(bnd_node,3),'LineWidth',2,'Color','b'); hold on;
nBnd_node = length(bnd_node);

%find bnd elements
ine=cell(np,1);
for j=1:ne
    for k=1:i34(j)
        nd=ele(j,k);
        ine{nd}=[ine{nd} j];
    end
end

bnd_ele=nan(nBnd_node-1,1); nBnd_ele=nBnd_node-1;
bnd_ele_n1=bnd_ele; bnd_ele_n2=bnd_ele;
for i=1:nBnd_node-1
    id1=bnd_node(i); id2=bnd_node(i+1);
    tmp = intersect(ine{id1},ine{id2});
    if (length(tmp)>1)
        bnd_ele(i) = tmp(1);
    else
        bnd_ele(i) = tmp;
    end
    bnd_ele_n1(i)=id1;
    bnd_ele_n2(i)=id2;
end
%sanity check
if sum(bnd_ele>0)~=nBnd_ele
    display('missing bnd ele');
    pause
end
%diagnostic plot
% scatter(node(bnd_ele_n1,2),node(bnd_ele_n1,3),5); hold on;
bnd_ele_x=nan(nBnd_ele,1); bnd_ele_y=bnd_ele_x;
for i=1:nBnd_ele
    el=bnd_ele(i);
    bnd_ele_x(i,1)=mean( node(ele(el,1:i34(el)) , 2 ) );
    bnd_ele_y(i,1)=mean( node(ele(el,1:i34(el)) , 3 ) );
end
%diagnostic plot
% scatter(bnd_ele_x,bnd_ele_y,5); hold on;

%-----(3) Read NWM segments info --------
%Note that we start from a sub-region ('shp_DB.shp') of the original
%NWM v1.2 shapefile to save time; the latter can also be used but it takes 
%much longer.
if exist([hdir 'X_NWM_SEG.shp'],'file')
    X_NWM_SEG=shaperead([hdir 'X_NWM_SEG.shp']);
else
    NWM_SEG=shaperead('shp_DB.shp');
    featureID=[NWM_SEG.featureID]';
    nseg=length(featureID);
    
    xx=node(bnd_node,2);
    yy=node(bnd_node,3);
    L2=[xx yy]';

    inter=cell(nseg,2);
    i_inter=zeros(nseg,1); n_inter=zeros(nseg,1);
    for i=1:nseg
        seg = NWM_SEG(i);
        L1 = [seg.X; seg.Y];
        P = InterX(L1,L2);
        if (~isempty(P))
            inter{i,1}=P(1,:);
            inter{i,2}=P(2,:);
            i_inter(i)=1;
            n_inter(i)=size(P,2);
        end
    end

    X_NWM_SEG=NWM_SEG(i_inter==1);     
    shapewrite(X_NWM_SEG,[hdir 'X_NWM_SEG.shp']);   
    
    %diagnostic plot
    P_all=[]; 
    for i=1:nseg
        if (~isempty(inter{i,1}))
            if (n_inter(i)>1) 
                display([num2str(n_inter(i)) ' intersections at: ' num2str(i)]);
            end
            for j=1:length(inter{i,1})
%                 scatter(inter{i,1},inter{i,2},5); hold on;
                P_all = [P_all  [inter{i,1};inter{i,2}]];
            end
        end
    end
end %if intersecting segs have already been found and saved in X_NWM_SEG.shp

n_X_Seg=length(X_NWM_SEG);
display(['Number of intersecting NWM segments: ' num2str(n_X_Seg)]);




%-----(4) find bnd elements and decide inflow/outflow--------
% Important: this assumes vertices on each segment are arranged from
% upstream to downstream, so the flow direction is determined
fname=[hdir 'NWM_to_SCHISM.mat'];
if exist(fname,'file')
    load(fname);
else
    P_inter = cell(n_X_Seg,5);
    for i=1:n_X_Seg
        display([num2str(i) ' of ' num2str(n_X_Seg)]);
        seg = X_NWM_SEG(i);
        for j=1:nBnd_ele
            L1 = [seg.X; seg.Y];
            L2=[node(bnd_ele_n1(j),2) node(bnd_ele_n2(j),2); 
                node(bnd_ele_n1(j),3) node(bnd_ele_n2(j),3)];
            P = InterX(L1,L2);

            %first screening, multiple intersections are possible for one
            %segment
            if (~isempty(P))
                %determine if intersection is already recorded in another element
                if ~isempty(P_inter{i,1})
                    for k=1:size(P,2)
                        tmp_dist = ( ( P_inter{i,1}-P(1,k) ).^2 + ( P_inter{i,2}-P(2,k) ).^2 ) .^0.5;
                        if min( tmp_dist )< 1e-10
                            if (size(P,2)==1)
                                P=[];
                            end 
                            P=[P(:,1:k-1) P(:,k+1:end)];
                        end
                    end
                    if size(P,2)==0 %no new intersections
                        break;
                    end
                end

                P_inter{i,1}=[P_inter{i,1},P(1,:)]; % intersection, x
                P_inter{i,2}=[P_inter{i,2},P(2,:)]; % intersection, y
                P_inter{i,3}=[P_inter{i,3} bnd_ele(j)]; %element
                P_inter{i,5}=[P_inter{i,5} size(P,2)]; %record number of intersections with this ele

		%Temporary (not really used): add a point in the ocean to close the polygon of the Delaware Bay
		%Use the intersecting element as the polygon in the future
                %in_poly=inpolygon(seg.X(1:end-1),seg.Y(1:end-1),[node(bnd_node,2); -73.685884], [node(bnd_node,3); 37.639312]);

                n_inout=0; %record in/out flow
                %second screening on each sub-segment (between two adjacent vertices
                %, only one intersection is possible
                for m=2:size(seg.X,2)-1
                    L10 = [seg.X(m-1:m); seg.Y(m-1:m)];
                    P0 = InterX(L10,L2);

                    if (size(P0,2)==1)

                        %nudge intersection toward seg node m by a little bit (1e-4 degree)
                        dp=([seg.X(m); seg.Y(m)] - P0) ;
                        dp_mag=sqrt(dp(1)^2+dp(2)^2);
                        P0_n=dp/dp_mag * 1e-4 + P0;
			%Use the intersecting element as the polygon in the future
                        if inpolygon(P0_n(1),P0_n(2),[node(bnd_node,2); -73.685884], [node(bnd_node,3); 37.639312])
                            n_inout=n_inout-1; %inflow
                        else
                            n_inout=n_inout+1; %outflow
                        end

                    elseif (size(P0,2)>1)
                        display('more than 1 intersection');
                        pause;  
                    end  
                end %for seg points                     
                P_inter{i,4}=[P_inter{i,4} n_inout]; %record in/outflow
            end %if (~isempty(P))

        end %for j=1:nBnd_ele
    end %for i=1:n_X_Seg
    save(fname,'P_inter');
end

% -----(5) read NWM outputs-----
TimeOri=datenum('1970-01-01');
RepoId=[X_NWM_SEG.featureID];

fname=[hdir 'NWM_flow.mat'];
if exist(fname,'file')
    load(fname);
    
    %tempo
%     RepoId=[2590277; 4784831]; %Trenton; Schu. R.

%     TM0=[];
%     SF0=[];
%     for  iday=TimeStart:1/24:TimeEnd
%         RepoName=[ REPODir datestr(iday,'yyyymmddHH') '00.CHRTOUT_DOMAIN1.comp'];
%         ncid0=netcdf.open(RepoName,'NC_NOWRITE');   
% 
%         %time
%         vid=netcdf.inqVarID(ncid0,'time');
%         timetmp=double(netcdf.getVar(ncid0,vid))/1440+TimeOri-TimeStart;
%         TM0=[TM0,timetmp];  
% 
%         %id list
%         vid=netcdf.inqVarID(ncid0,'feature_id');
%         id=double(netcdf.getVar(ncid0,vid));
%         [tmp  I2]=ismember(RepoId,id);
%          
%         %stream
%         vid=netcdf.inqVarID(ncid0,'streamflow');
%         sftmp=double(netcdf.getVar(ncid0,vid))*0.01;
%         SF0=[SF0, sftmp(I2)];    
%     end
else
    
    TM=[];
    QL=[];
    SF=[];
    for  iday=TimeStart:1/24:TimeEnd
        RepoName=[ REPODir datestr(iday,'yyyymmddHH') '00.CHRTOUT_DOMAIN1.comp']
        ncid0=netcdf.open(RepoName,'NC_NOWRITE');   

        %time
        vid=netcdf.inqVarID(ncid0,'time');
        timetmp=double(netcdf.getVar(ncid0,vid))/1440+TimeOri-TimeStart;
        TM=[TM,timetmp];  

        %id list
        vid=netcdf.inqVarID(ncid0,'feature_id');
        id=double(netcdf.getVar(ncid0,vid));
        [tmp  I2]=ismember(RepoId,id);
        nMissing=sum(I2==0);
        if iday==TimeStart && nMissing>0
            display(['warning: ' num2str(nMissing) ' intersecting segs are not found in nc data, featureId(s):']);
            X_NWM_SEG(find(tmp==0)).featureID
            display('press any key to continue or ctrl-c to abort ...');
            pause
            RepoId=RepoId(I2~=0); P_inter=P_inter(I2~=0,:); I2=I2(I2~=0);
            X_NWM_SEG=X_NWM_SEG((tmp==1));
        end

        %lateral
        vid=netcdf.inqVarID(ncid0,'q_lateral');
        qltmp=double(netcdf.getVar(ncid0,vid))*0.1;
        QL= [QL, qltmp(I2)];
        %stream
        vid=netcdf.inqVarID(ncid0,'streamflow');
        sftmp=double(netcdf.getVar(ncid0,vid))*0.01;
        SF=[SF, sftmp(I2)];   
        
        netcdf.close(ncid0);
    end
    save(fname,'SF','QL','X_NWM_SEG','RepoId','P_inter');
end

%!!!!!!skip upstream 
% for i=1:length(X_NWM_SEG)
%     if X_NWM_SEG(i).lat>40.2483
%         QL(i,:)=0.0; SF(i,:)=0.0;
%     end
% end

% -----(6) make SCHISM input files-----
% T is specified as -9999,
% which is an identifier in SCHISM meaning T is not specified in in/outflow
% and takes the ambient value.


source_ele=[]; sink_ele=[];
source_flow=[]; sink_flow=[];
nsource=0; nsink=0;
for i=1:size(P_inter,1)
    %tempo
    if ismember(X_NWM_SEG(i).featureID,skip_segs)
        continue;
    end
     
    for j=1:length(P_inter{i,4})
        if (P_inter{i,4}(j)==-1)
            nsource=nsource+1;
            source_ele=[source_ele P_inter{i,3}(j)];
            source_flow=[source_flow; SF(i,:)];
        elseif (P_inter{i,4}(j)==1)
            nsink=nsink+1;
            sink_ele=[sink_ele P_inter{i,3}(j)];
            sink_flow=[sink_flow; SF(i,:)];
        end
    end
end

unique_source_ele=[]; unique_source_flow=[];
for i=1:nsource
    [ie id]=ismember(source_ele(i),unique_source_ele);
    if (id==0)
        unique_source_ele=[unique_source_ele source_ele(i)];
        unique_source_flow=[unique_source_flow; source_flow(i,:)];
    else
        unique_source_flow(id,:)=unique_source_flow(id,:)+source_flow(i,:);
    end
end

unique_sink_ele=[]; unique_sink_flow=[];
for i=1:nsink
    [ie id]=ismember(sink_ele(i),unique_sink_ele);
    if (id==0)
        unique_sink_ele=[unique_sink_ele sink_ele(i)];
        unique_sink_flow=[unique_sink_flow; sink_flow(i,:)];
    else
        unique_sink_flow(id,:)=unique_sink_flow(id,:)+sink_flow(i,:);
    end
end
    

fid=fopen([hdir 'source_sink.in'],'wt');
fprintf(fid,'%d\n',length(unique_source_ele)+n_plus); %plus Susq. R
for i=1:length(unique_source_ele)
    fprintf(fid,'%d\n',unique_source_ele(i));
end
fprintf(fid,'%d\n',ie_plus); %Susq.

fprintf(fid,'\n');
fprintf(fid,'%d\n',length(unique_sink_ele));
for i=1:length(unique_sink_ele)
    fprintf(fid,'%d\n',unique_sink_ele(i));
end
fclose(fid);

time_stamp=(0:3600:(TimeEnd-TimeStart)*86400);
susq_r=susq_r_flow*ones(size(time_stamp));
% trenton=load('F:\Projects\Delaware\NWM_Coupling\usgs_01463500_flow.DATA');
% trenton=trenton(1:4:1201*4)';
dlmwrite([hdir 'vsource.th'],[time_stamp; unique_source_flow; susq_r]','delimiter',' ','precision',15);

salt=zeros(size(unique_source_flow,1)+n_plus,size(unique_source_flow,2));
temp=-9999*ones(size(unique_source_flow,1)+n_plus,size(unique_source_flow,2));
dlmwrite([hdir 'msource.th'],[time_stamp; temp; salt]','delimiter',' ','precision',15);

dlmwrite([hdir 'vsink.th'],[time_stamp; -unique_sink_flow]','delimiter',' ','precision',15);

viz_source;






















