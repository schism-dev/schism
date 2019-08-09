function path = partpath(pthfile)

[fid,error] = fopen(pthfile);
if fid<0
    disp(error)
end
[a,count] = fscanf(fid,'%s',1);
[a,count] = fscanf(fid,'%i',1);

step = 0;
while count>0
    step = step+1;
    [a,count] = fscanf(fid,'%f %i',2);
    %if count==0 | step==4
    if count==0
        path.p = p;
        path.q = q;
        path.x = x;
        path.y = y;
        path.z = z;
        path.wx=wx;
        path.wy=wy;
        path.el=el;
        path.pnum=pnum;
        path.time = time;
        return
    end
    time(step) = a(1);
    pnum(step) = a(2);

    for i = 1:pnum(step)
        [data,count] = fscanf(fid,'%i %i %f %f %f %f',6);
        if count~=6
            break
        end
        p(i,step) = data(2);
        q(i,step) = data(3);
        x(i,step) = data(4);
        y(i,step) = data(5);
        z(i,step) = data(6);
    end

     [a,count] = fscanf(fid,'%f %f %f',3);
     wx(step)=a(1);
     wy(step)=a(2);
     el(step)=a(3);
end

