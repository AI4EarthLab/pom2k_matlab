depth();
Problem_name=input('Please input the name of problem (seamount, box or file2ic):','s');
while (~strcmp(Problem_name,'seamount') &&~strcmp(Problem_name,'box') &&~strcmp(Problem_name,'file2ic'))
    Problem_name=input('Input error!Try again:','s');
end
if(strcmp(Problem_name,'seamount'))
    [tb,sb,tclim,sclim,wusurf,wvsurf,wtsurf,wssurf,uvel,vvel,swrad] = seamount(e_atmos);
elseif(strcmp(Problem_name,'box'))
    [tb,sb,tclim,sclim,wusurf,wvsurf,wtsurf,wssurf,uvel,vvel,swrad,vfluxf,e_atmos] = box(e_atmos);
else
    file2ic_name=input('Please input the name of file2ic (For example:XX.dat):','s');
    [tb,sb,tclim,sclim,wusurf,wvsurf,wtsurf,wssurf,uvel,vvel,swrad,vfluxf,e_atmos] = file2ic(file2ic_name);
end
