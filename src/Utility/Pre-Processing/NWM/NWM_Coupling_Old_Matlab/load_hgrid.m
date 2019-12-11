function [ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(varargin)

    if nargin==0
        return;
    else
        hdir=varargin{1,1};
        fname1=varargin{1,2};
        fname2=varargin{1,3};
        ibnd=varargin{1,4};
    end
    
    
    [ne, np]=textread([hdir '/' fname1], '%d %d', 1, 'headerlines', 1);
    i34=zeros(ne,1); ele=nan.*ones(ne,4);
    
    [node(:,1),node(:,2),node(:,3),node(:,4)]=textread([hdir '/' fname2], '%d %f %f %f',np, 'headerlines', 2);
    
    fid = fopen([hdir  '/' fname1]);

    for i=1:2+np
        tline = fgets(fid);
    end
    for i=1:ne
        tline = strtrim(fgets(fid));
        tmp=strsplit(tline);
        i34(i)=str2num(tmp{2});
        for j=1:i34(i)
            ele(i,j)=str2num(tmp{2+j});
        end
    end

    fclose(fid);
    
    if ibnd==0
        bndnode=[];
        return;
    end
    
    %bnd
    bndnode=[];
    
    ptr=np+ne+2;
    [nope]=textread([hdir '/' fname1], '%d',1, 'headerlines', ptr);  ptr=ptr+1;  
    [dummy]=textread([hdir  '/' fname1], '%d',1, 'headerlines', ptr);  ptr=ptr+1;
    open_bnds=cell(nope,1);
    for i=1:nope
        [nob(i)]=textread([hdir '/' fname1], '%d',1, 'headerlines', ptr); ptr=ptr+1;
        [tmp]=textread([hdir '/' fname1], '%d',nob(i), 'headerlines', ptr); ptr=ptr+nob(i);
        if i==1 
            bndnode=[bndnode; tmp];
        end
        open_bnds{i}=tmp;
    end
    
    ocean0=[];
    [nlbe]=textread([hdir  '/' fname1], '%d',1, 'headerlines', ptr);  ptr=ptr+1;  
    [dummy]=textread([hdir  '/' fname1], '%d',1, 'headerlines', ptr);  ptr=ptr+1;   
    land_bnds=cell(nlbe,1); ilb_island=zeros(nlbe,1);
    for i=1:nlbe
        [nlb(i) ilb_island(i)]=textread([hdir '/' fname1], '%d %d',1, 'headerlines', ptr); ptr=ptr+1;
        [tmp]=textread([hdir  '/' fname1], '%d',nlb(i), 'headerlines', ptr); ptr=ptr+nlb(i);    
        if (ilb_island(i)==1 && isempty(ocean0))
            ocean0=bndnode;
        end
        bndnode=[bndnode; tmp];
        land_bnds{i}=tmp;
    end
        
    bndnode(:,2:4)=node(bndnode(:,1),2:4);
    ocean=node(ocean0(:,1),2:3);
    

    if ibnd==2
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
    end
end
    