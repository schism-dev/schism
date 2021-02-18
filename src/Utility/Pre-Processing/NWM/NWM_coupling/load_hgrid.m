function [ne,np,node,ele,i34,bndnode,open_bnds,land_bnds,ilb_island]=load_hgrid(varargin)
    
    ele=[]; i34=[]; bndnode=[]; open_bnds=[]; land_bnds=[]; ilb_island=[];

    if nargin==0
        return;
    else
        hdir=varargin{1,1};
        fname1=varargin{1,2};
        fname2=varargin{1,3};
        ibnd=varargin{1,4};
    end
    
    fid=fopen([hdir '/' fname1]);
    tmp=textscan(fid, '%d%d', 1, 'headerlines', 1);
    ne=tmp{1}(1); np=tmp{2}(1);
    i34=zeros(ne,1); ele=nan(ne,4); node=nan(np,4);
    
    tmp=textscan(fid, '%d %f %f %f',np);
    node(:,2:4)=cell2mat(tmp(2:4));
    
    if ibnd==-1 %not even reading the elements
        fclose(fid);
        return
    end
    
    tmp=textscan(fid, '%d %d %d %d %d %d',ne);
    ele(:,1:4)=cell2mat(tmp(3:end)); ele(ele==0)=nan;
    i34=cell2mat(tmp(2));
        
    if ibnd==0
        bndnode=[];
        fclose(fid);
        return;
    end
    
    %bnd
    bndnode=[];
    
    tmp=textscan(fid,'%d',1); dummy = fgetl(fid); %get the rest of the line
    nope=tmp{1};
    
    tmp=textscan(fid,'%d',1); dummy = fgetl(fid); %get the rest of the line
    n_op_nodes=tmp{1};
    
    open_bnds=cell(nope,1);
    for i=1:nope
        tmp=textscan(fid,'%d',1); dummy = fgetl(fid); %get the rest of the line
        nob(i)=tmp{1};        
        
        tmp=textscan(fid,'%d',nob(i));
        tmp=tmp{1};
        if i==1 
            bndnode=[bndnode; tmp];
        end
        open_bnds{i}=tmp;
    end
        
    tmp=textscan(fid,'%d',1); dummy = fgetl(fid); %get the rest of the line
    nlbe=tmp{1};    
    tmp=textscan(fid,'%d',1); dummy = fgetl(fid); %get the rest of the line
    n_lb_nodes=tmp{1};
    
    land_bnds=cell(nlbe,1); ilb_island=zeros(nlbe,1);
    for i=1:nlbe
        tmp=textscan(fid,'%d',1); dummy = fgetl(fid); %get the rest of the line
        nlb(i)=tmp{1}; ilb_island(i)=tmp{1};
        tmp=textscan(fid,'%d',nlb(i));
        tmp=tmp{1};
        
        bndnode=[bndnode; tmp];
        land_bnds{i}=tmp;
    end
        
    bndnode(:,2:4)=node(bndnode(:,1),2:4);    

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
    
    fclose(fid);
end
    