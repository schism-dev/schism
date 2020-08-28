clear;
%------------------usage-----------------------
%  Important: only works for T=-9999 and S=0
%  Make sure wDir has hgrid.ll and two sets of source/sink:
%   source_sink.in.[12]
%   vsource.th.[12]
%   msource.th.[12]
%   optional: vsink.th.[12]
%  The outputs are source/sink files with regular names
wDir='/sciclone/schism10/feiye/work/Gulf_Stream/RUN36/';
%----------------------------------------------

[ne,np]=textread([wDir 'hgrid.ll'], '%d%d',1, 'headerlines', 1);

soid=cell(2,1); siid=soid; vso=soid; vsi=soid;

display('reading existing source/sinks');
for i=1:2
    
    %read source ele
    [ne_source(i)]=textread([wDir 'source_sink.in.' num2str(i)], '%d', 1, 'headerlines', 0);
    [soid{i}]=textread([wDir 'source_sink.in.' num2str(i)], '%d', ne_source(i), 'headerlines', 1);
    vso{i}=load([wDir 'vsource.th.' num2str(i)]);

    %read sink ele if any
    [ne_sink(i)]=textread([wDir 'source_sink.in.' num2str(i)], '%d', 1, 'headerlines', 1+ne_source(i)+1);
    if ne_sink(i)>0
        [siid{i}]=textread([wDir 'source_sink.in.' num2str(i)], '%d', ne_sink, 'headerlines',1+ne_source(i)+1+1);
        vsi{i}=load([wDir 'vsink.th.' num2str(i)]);
    end

    dt(i)=vso{i}(2,1)-vso{i}(1,1);
    last_t(i)=vso{i}(end,1);
end


last_t0=min(last_t);
dt0=min(dt);

%common ie
ioSS=zeros(ne,1); iiSS=ioSS;
for i=1:2
    ioSS(soid{i})=1;
    if (ne_sink(i)>0)
        iiSS(siid{i})=1;
    end
end

iso=find(ioSS==1);
isi=find(iiSS==1);

time_stamp=[0:dt0:last_t0]'; nt0=length(time_stamp);

%combine sources
display('combining sources');
vs_interp=cell(2,1);
for i=1:2
    vs_interp{i}=nan(nt0,length(iso));
    tmp=interp1(vso{i}(:,1), vso{i}(:,2:end),time_stamp);
    
    %find index in the new source_sink.in
    [tf,loc]=ismember(soid{i},iso);
    if (min(tf)==0)
        display('error: index not found');
        return;
    end
    vs_interp{i}(:,loc)=tmp;        
end
    
    
II=find(isnan(vs_interp{1}(1,:))&isnan(vs_interp{2}(1,:)));
if ~isempty(II)
    display('error: nan found in vsource_interp');
    return;
end

for i=1:2
    vs_interp{i}(isnan(vs_interp{i}))=0;
end

vs_combined=vs_interp{1}+vs_interp{2};
dlmwrite([wDir 'vsource.th'],[time_stamp vs_combined],'precision',10,'delimiter',' ');
clearvars vs_combined vs_interp;

if (max(ne_sink)>0) 
  display('combining sinks');
    %combine sinks
    vs_interp=cell(2,1);
    for i=1:2
        if ne_sink(i)>0 
            vs_interp{i}=nan(nt0,length(isi));
            tmp=interp1(vsi{i}(:,1), vsi{i}(:,2:end),time_stamp);

            %find index in the new source_sink.in
            [tf,loc]=ismember(siid{i},isi);
            if (min(tf)==0)
                display('error: index not found');
                return;
            end
            vs_interp{i}(:,loc)=tmp;        
        end
    end


    if (ne_sink(1)>0 && ne_sink(2)==0)
        II=find(isnan(vs_interp{1}(1,:)));
    elseif (ne_sink(1)==0 && ne_sink(2)>0)
        II=find(isnan(vs_interp{2}(1,:)));
    else
        II=find(isnan(vs_interp{1}(1,:))&isnan(vs_interp{2}(1,:)));
    end
    if ~isempty(II)
        display('error: nan found in vsink_interp');
        return;
    end

    vs_combined=zeros(nt0,length(isi));
    for i=1:2
        if ne_sink(i)>0
            vs_interp{i}(isnan(vs_interp{i}))=0;
            vs_combined=vs_combined+vs_interp{i};
        end
    end

    dlmwrite([wDir 'vsink.th'],[time_stamp vs_combined],'precision',10,'delimiter',' ');
end %if any sink        

%write combined source_sink.in
display('writing other files');
fid=fopen([wDir '/source_sink.in'],'wt');
fprintf(fid,'%d\n',length(iso)); %number of sources
for i=1:length(iso)
    fprintf(fid,'%d\n',iso(i));
end
fprintf(fid,'\n');
fprintf(fid,'%d\n',length(isi)); %number of sinks
for i=1:length(isi)
    fprintf(fid,'%d\n',isi(i));
end
fclose(fid);

%write new msource
total_nt=length(time_stamp);
msource=zeros(total_nt,length(iso)*2+1);
msource(:,1)=time_stamp;
msource(:,2:length(iso)+1)=-9999;
dlmwrite([wDir 'msource.th'],msource,'precision',10,'delimiter',' ');
            
    
