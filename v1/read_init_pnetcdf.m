function [east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,tb,sb,tclim, ...
       sclim,rot,uvel,vvel,vfluxf,wusurf,wvsurf,e_atmos,ub,vb,uab,vab, ...
        elb,etb,dt,uabw,uabe,vabs,vabn,eln,els,ele,elw,ssurf,tsurf,tbe,tbw,tbn,tbs, ...
                           sbe,sbw,sbn,sbs] = read_init_pnetcdf
global z zz dx dy cor h fsm dum dvm im kb jm z_3d zz_3d dz dzz art aru arv ...
dx_3d dy_3d dz_3d dzz_3d cor_3d h_3d art_3d aru_3d arv_3d fsm_3d dum_3d dvm_3d rfe rfw rfn rfs problem in_path;
fnc=[in_path,problem,'.nc'];
% fnc='seamount.nc';
% fnc=input('Enter the path and name of the nc file:','s');
% while (~exist(fnc))
%     disp('The file is not exist!');
%     fnc = input('Please enter again:','s');
% end
% Get the values

z       = ncread(fnc,'z')        ;
zz      = ncread(fnc,'zz')       ;
dz      = ncread(fnc,'dz')       ;
dzz     = ncread(fnc,'dzz')      ;
dx      = ncread(fnc,'dx')       ;
dy      = ncread(fnc,'dy')       ;
cor     = ncread(fnc,'cor')      ;
h       = ncread(fnc,'h')        ;
fsm     = ncread(fnc,'fsm')      ;
dum     = ncread(fnc,'dum')      ;
dvm     = ncread(fnc,'dvm')      ;
art     = ncread(fnc,'art')      ;
aru     = ncread(fnc,'aru')      ;
arv     = ncread(fnc,'arv')      ;
rfe     = ncread(fnc,'rfe')      ;
rfw     = ncread(fnc,'rfw')      ;
rfn     = ncread(fnc,'rfn')      ;
rfs     = ncread(fnc,'rfs')      ;
east_e  = ncread(fnc,'east_e')   ;
north_e = ncread(fnc,'north_e')  ;
east_c  = ncread(fnc,'east_c')   ;
north_c = ncread(fnc,'north_c')  ;
east_u  = ncread(fnc,'east_u')   ;
north_u = ncread(fnc,'north_u')  ;
east_v  = ncread(fnc,'east_v')   ;
north_v = ncread(fnc,'north_v')  ;
tb      = ncread(fnc,'tb')       ;
sb      = ncread(fnc,'sb')       ;
tclim   = ncread(fnc,'tclim')    ;
sclim   = ncread(fnc,'sclim')    ;
rot     = ncread(fnc,'rot')      ;
uvel    = ncread(fnc,'uvel')     ;
vvel    = ncread(fnc,'vvel')     ;
vfluxf  = ncread(fnc,'vfluxf')   ;
wusurf  = ncread(fnc,'wusurf')   ;
wvsurf  = ncread(fnc,'wvsurf')   ;
e_atmos = ncread(fnc,'e_atmos')  ;
ub      = ncread(fnc,'ub')       ;
vb      = ncread(fnc,'vb')       ;
uab     = ncread(fnc,'uab')      ;
vab     = ncread(fnc,'vab')      ;
elb     = ncread(fnc,'elb')      ;
etb     = ncread(fnc,'etb')      ;
dt      = ncread(fnc,'dt')       ;
uabw    = ncread(fnc,'uabw')     ;
uabe    = ncread(fnc,'uabe')     ;
vabs    = ncread(fnc,'vabs')     ;
vabn    = ncread(fnc,'vabn')     ;
els     = ncread(fnc,'els')      ;
eln     = ncread(fnc,'eln')      ;
ele     = ncread(fnc,'ele')      ;
elw     = ncread(fnc,'elw')      ;
ssurf   = ncread(fnc,'ssurf')    ;
tsurf   = ncread(fnc,'tsurf')    ;
tbe     = ncread(fnc,'tbe')      ;
sbe     = ncread(fnc,'sbe')      ;
sbw     = ncread(fnc,'sbw')      ;
tbw     = ncread(fnc,'tbw')      ;
tbn     = ncread(fnc,'tbn')      ;
tbs     = ncread(fnc,'tbs')      ;
sbn     = ncread(fnc,'sbn')      ;
sbs     = ncread(fnc,'sbs')      ;
 for k=1:kb
  z_3d(:,:,k) = repmat(z(k),im,jm);   zz_3d(:,:,k) = repmat(zz(k),im,jm);
  dz_3d(:,:,k)= repmat(dz(k),im,jm);  dzz_3d(:,:,k) =   repmat(dzz(k),im,jm);
 end
 
    dx_3d=repmat(dx,1,1,kb);    dy_3d=repmat(dy,1,1,kb);
    cor_3d=repmat(cor,1,1,kb);  h_3d=repmat(h,1,1,kb);
    art_3d=repmat(art,1,1,kb);  aru_3d=repmat(aru,1,1,kb);
    arv_3d=repmat(arv,1,1,kb);  fsm_3d=repmat(fsm,1,1,kb);
    dum_3d=repmat(dum,1,1,kb);  dvm_3d=repmat(dvm,1,1,kb);  
    
end