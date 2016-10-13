function write_to_file(im,jm,kb,time,wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,...
    utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,...
    km,kh,kq,l,q2,q2b,aam,q2l,q2lb)



fid = fopen('data.out','wt');
fprintf(fid,'%f',time);
for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',wubot(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',wvbot(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',aam2d(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',ua(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',uab(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',va(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',vab(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',el(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',elb(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',et(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',etb(i,j));
    end
    fprintf(fid,'\n');
end


for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',egb(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',utb(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',vtb(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',u(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf('\n');
end
for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',ub(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf('\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',w(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end


for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',v(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',vb(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',t(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',tb(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',s(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',sb(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',rho(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end


for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',adx2d(i,j));
    end
    fprintf(fid,'\n');
end


for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',ady2d(i,j));
    end
    fprintf(fid,'\n');
end


for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',advua(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        fprintf(fid,'%f ',advva(i,j));
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',km(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',kh(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',kq(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',l(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',q2(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',q2b(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',aam(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',q2l(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

for i=1:im
    for j=1:jm
        for k=1:kb
            fprintf(fid,'%f ',q2lb(i,j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

fclose(fid);
return
end
