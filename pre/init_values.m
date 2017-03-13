depth();
Problem_name=input('Please input the name of problem (seamount, box or file2ic):','s');
while (~strcmp(Problem_name,'seamount') &&~strcmp(Problem_name,'box') &&~strcmp(Problem_name,'file2ic'))
    Problem_name=input('Input error!Try again:','s');
end
if(strcmp(Problem_name,'seamount'))
    [tb,sb,tclim,sclim,ub,uab,elb,etb,dt,aam2d,rho,rmean,wusurf,wvsurf,dt_3d] = seamount(e_atmos,aam);
elseif(strcmp(Problem_name,'box'))
    [tb,sb,tclim,sclim,ub,uab,elb,etb,aam2d,dt,dt_3d,rho,rmean,tatm,satm,vfluxf,e_atmos] = box(e_atmos,aam);
else
    file2ic_name=input('Please input the name of file2ic (For example:XX.dat):','s');
    [tb,sb,tclim,sclim,ub,vb,elb,etb,dt,rho,rmean,ssurf,tsurf,dt_3d] = file2ic(file2ic_name);
end
