function [dx,dy,z,zz,cor,rot,h,fsm,dum,dvm,t,s,tclim,sclim,wusurf,wvsurf,wtsurf,  ...
            wssurf,swrad,vflux,e_atmos]=read_basic_variables(problem)
   fnc      =['basic_file/',problem,'_basic.nc']    ;
   dx       = ncread(fnc,'dx')        ;
   dy       = ncread(fnc,'dy')        ;
   z        = ncread(fnc,'z')         ;
   zz       = ncread(fnc,'zz')        ;
   cor      = ncread(fnc,'cor')       ;
   rot      = ncread(fnc,'rot')       ;
   h        = ncread(fnc,'h')         ;
   fsm      = ncread(fnc,'fsm')       ;
   dum      = ncread(fnc,'dum')       ;
   dvm      = ncread(fnc,'dvm')       ;
   t        = ncread(fnc,'t')         ;
   s        = ncread(fnc,'s')         ;
   tclim    = ncread(fnc,'tclim')     ;
   sclim    = ncread(fnc,'sclim')     ;
   wusurf   = ncread(fnc,'wusurf')    ;
   wvsurf   = ncread(fnc,'wvsurf')    ;
   wtsurf   = ncread(fnc,'wtsurf')    ;
   wssurf   = ncread(fnc,'wssurf')    ;
   swrad    = ncread(fnc,'swrad')     ;
   vflux    = ncread(fnc,'vflux')     ;
   e_atmos  = ncread(fnc,'e_atmos')   ;
end