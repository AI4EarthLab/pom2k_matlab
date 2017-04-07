if (strcmp(Problem_name,'seamount'))
    ncid      =    netcdf.create('seamount_basic.nc','NC_SHARE');
elseif (strcmp(Problem_name,'box'))
    ncid      =    netcdf.create('box_basic.nc','NC_SHARE');
else
    ncid      =    netcdf.create('file2ic_basic.nc','NC_SHARE');   
end
imx       =    netcdf.defDim(ncid,'im',im);     jmy       =    netcdf.defDim(ncid,'jm',jm);
kbz       =    netcdf.defDim(ncid,'kb',kb);     dconst    =    netcdf.defDim(ncid,'dc',1);


zid       =    netcdf.defVar(ncid,'z','double',[dconst kbz]);
zzid      =    netcdf.defVar(ncid,'zz','double',[dconst kbz]);
dxid      =    netcdf.defVar(ncid,'dx','double',[imx jmy]);
dyid      =    netcdf.defVar(ncid,'dy','double',[imx jmy]);
corid     =    netcdf.defVar(ncid,'cor','double',[imx jmy]);
hid       =    netcdf.defVar(ncid,'h','double',[imx jmy]);
fsmid     =    netcdf.defVar(ncid,'fsm','double',[imx jmy]);
dumid     =    netcdf.defVar(ncid,'dum','double',[imx jmy]);
dvmid     =    netcdf.defVar(ncid,'dvm','double',[imx jmy]);
rotid     =    netcdf.defVar(ncid,'rot','double',[imx jmy]);
wusurfid  =    netcdf.defVar(ncid,'wusurf','double',[imx jmy]);
wvsurfid  =    netcdf.defVar(ncid,'wvsurf','double',[imx jmy]);
wtsurfid  =    netcdf.defVar(ncid,'wtsurf','double',[imx jmy]);
wssurfid  =    netcdf.defVar(ncid,'wssurf','double',[imx jmy]);
e_atmosid =    netcdf.defVar(ncid,'e_atmos','double',[imx jmy]);
swradid   =    netcdf.defVar(ncid,'swrad','double',[imx jmy]);
vfluxfid  =    netcdf.defVar(ncid,'vflux','double',[imx jmy]);
tbid     =    netcdf.defVar(ncid,'t','double',[imx jmy kbz]);
sbid     =    netcdf.defVar(ncid,'s','double',[imx jmy kbz]);
tclimid  =    netcdf.defVar(ncid,'tclim','double',[imx jmy kbz]);
sclimid  =    netcdf.defVar(ncid,'sclim','double',[imx jmy kbz]);

netcdf.endDef(ncid)

netcdf.putVar(ncid,zid,z);              netcdf.putVar(ncid,zzid,zz);
netcdf.putVar(ncid,dxid,dx);            netcdf.putVar(ncid,dyid,dy);
netcdf.putVar(ncid,corid,cor);          netcdf.putVar(ncid,hid,h);      
netcdf.putVar(ncid,fsmid,fsm);          netcdf.putVar(ncid,rotid,rot);
netcdf.putVar(ncid,dumid,dum);          netcdf.putVar(ncid,dvmid,dvm);
netcdf.putVar(ncid,wusurfid,wusurf);      netcdf.putVar(ncid,wvsurfid,wvsurf);    
netcdf.putVar(ncid,wtsurfid,wtsurf);      netcdf.putVar(ncid,wssurfid,wssurf);       
netcdf.putVar(ncid,tbid,tb);              netcdf.putVar(ncid,sbid,sb);
netcdf.putVar(ncid,tclimid,tclim);        netcdf.putVar(ncid,sclimid,sclim);
netcdf.putVar(ncid,e_atmosid,e_atmos);    netcdf.putVar(ncid,swradid,swrad);
netcdf.putVar(ncid,vfluxfid,vfluxf);
netcdf.close(ncid);