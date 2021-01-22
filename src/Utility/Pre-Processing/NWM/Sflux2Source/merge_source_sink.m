clear;
%------------------usage-----------------------
%  Important: only works for T=-9999 and S=0
%  Make sure wDir has hgrid.ll and two sets of source/sink:
%   source_sink.in.[12]
%   vsource.th.[12]
%   optional: vsink.th.[12]
%  The outputs are source/sink files with regular names:
%   source_sink.in
%   vsource.th
%   msource.th (T=-9999; S=0)
%   vsink.th

%------------------inputs-----------------------
wDir='./';
nday=1; %Days needed; recommended: 1 day more than the rnday in param.nml
%----------------------------------------------



[ne,np]=textread([wDir 'hgrid.ll'], '%d%d',1, 'headerlines', 1);

soid=cell(2,1); siid=soid; vso=soid; vsi=soid;
dt=nan(2,2);

for i=1:2    
    %read source ele
    [ne_source(i)]=textread([wDir 'source_sink.in.' num2str(i)], '%d', 1, 'headerlines', 0);
    [soid{i}]=textread([wDir 'source_sink.in.' num2str(i)], '%d', ne_source(i), 'headerlines', 1);
    fid=fopen([wDir 'vsource.th.' num2str(i)]);
    tmp0=textscan(fid, '%f', 1+ne_source(i));
    tmp1=textscan(fid, '%f', 1+ne_source(i));
    fclose(fid);
    dt(i,1)=tmp1{1}(1)-tmp0{1}(1); %in seconds

    %read sink ele if any
    [ne_sink(i)]=textread([wDir 'source_sink.in.' num2str(i)], '%d', 1, 'headerlines', 1+ne_source(i)+1);
    if ne_sink(i)>0
        [siid{i}]=textread([wDir 'source_sink.in.' num2str(i)], '%d', ne_sink(i), 'headerlines',1+ne_source(i)+1+1);
        fid=fopen([wDir 'vsource.th.' num2str(i)]);
        tmp0=textscan(fid, '%f', 1+ne_source(i));
        tmp1=textscan(fid, '%f', 1+ne_source(i));
        fclose(fid);
        dt(i,2)=tmp1{1}(1)-tmp0{1}(1); %in seconds
    end
end



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

fout=fopen([wDir '/source_sink.in'],'wt');
fprintf(fid,'%d\n',length(iso));
fprintf(fid,'%d\n',iso);
fprintf(fid,'\n');
fprintf(fid,'%d\n',length(isi));
fprintf(fid,'%d\n',isi);
fclose(fout);

%---------------------------------------------------------
%----------merge sources and sinks----------------------
%---------------------------------------------------------

fnames={'vsource.th','vsink.th'};
for iter=1:2
    
    dt0=min(dt(:,iter)); i_self_dt=dt0==dt(:,iter);
    time_stamp=[0:dt0:nday*86400]'; 
    nt0=length(time_stamp);
    
    if iter==1
        %map from [12] to combined
        ib=cell(2,1);
        [~, ib{1}]=ismember(soid{1},iso);
        [~, ib{2}]=ismember(soid{2},iso);
        ne_ss=ne_source;
        nrec=length(iso);
    else
        if ne_sink(1)==0 || ne_sink(2)==0
            if ne_sink(1)==0 && ne_sink(2)==0
              return
            end
            if ne_sink(1)>0
                copyfile 'vsink.th.1' 'vsink.th';
            elseif ne_sink(2)>0
                copyfile 'vsink.th.2' 'vsink.th';
            end
            return
        end
        %map from [12] to combined
        ib=cell(2,1);
        [~, ib{1}]=ismember(siid{1},isi);
        [~, ib{2}]=ismember(siid{2},isi);
        ne_ss=ne_sink;
        nrec=length(isi);
    end

    %one line of outputs
    fid=cell(2,1);
    fid{1}=fopen([wDir '/' fnames{iter} '.1']);
    fid{2}=fopen([wDir '/' fnames{iter} '.2']);
    delete([wDir '/' fnames{iter}]);

    %read two lines from each file to start the loop
    this_time = 0;
    tmp=cell(2,2); t=zeros(2,2);
    for i=1:2
        tmp(i,1)=textscan(fid{i}, '%f', 1+ne_ss(i)); 
        t(i,1)=tmp{i,1}(1);
        tmp(i,2)=textscan(fid{i}, '%f', 1+ne_ss(i)); 
        t(i,2)=tmp{i,2}(1);
    end
    while this_time < time_stamp(end)
        out = zeros(1,nrec);
        for i=1:2
            if t(i,2)< this_time+dt0
                tmp(i,1)=tmp(i,2);
                tmp(i,2)=textscan(fid{i}, '%f', 1+ne_ss(i));
                t(i,1)=t(i,2);
                t(i,2)=tmp{i,2}(1);
            end
            if i_self_dt(i)
                out(ib{i})=out(ib{i})+tmp{1,1}(2:end)';
            else
                tmp_val = cell2mat(tmp(i,1:2));
                tmp0 = interp1(t(i,1:2)',tmp_val(2:end,:)',this_time);
                out(ib{i})=out(ib{i})+tmp0;
            end
        end
        dlmwrite([wDir '/' fnames{iter}], [this_time out],'precision',10,'delimiter',' ','-append');
        this_time=this_time+dt0;
        waitbar(this_time / time_stamp(end))
    end
end


%---------------------------------------------------------
%           output dummy msource.th
%---------------------------------------------------------

msource=zeros(2,length(iso)*2+1);
msource(:,1)=[0 time_stamp(end)*1.1];
msource(:,2:length(iso)+1)=-9999;
dlmwrite([wDir 'msource.th'],msource,'precision',15,'delimiter',' ');


