[east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,tb,sb,tclim, ...
       sclim,rot,uvel,vvel,vfluxf,wusurf,wvsurf,e_atmos,ub,vb,uab,vab, ...
        elb,etb,dt,uabw,uabe,vabs,vabn,eln,els,ele,elw,ssurf,tsurf,tbe,tbw,tbn,tbs, ...
                           sbe,sbw,sbn,sbs] = read_init_pnetcdf;
% %Initialize density
[rho]=dens(sb,tb,h_3d,fsm_3d);
if fclim_flag
    [rmean]=dens(sclim,tclim,h_3d,fsm_3d);
else
    rmean(:,:,1:kbm1) = rho(:,:,1:kbm1);
end