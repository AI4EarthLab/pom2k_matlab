clear all;
init_global();
init_constants();

initialize_arrays();
read_input();

OP=init_operator();
gs=init_grid(gridtype);

init_fields();

    [ua,va,el,et,etf,d,dt,w,d_3d,dt_3d,l,q2b,q2lb,kh,km,kq,aam,q2,q2l,t,s,u,v,rho,drhox,drhoy,drx2d,dry2d] ...
          = update_initial(uab,vab,elb,etb,h,q2b,small,aam,aam_init,w,vfluxf,tb,sb,ub,vb,rmean,ramp);
    
    [cbc] = bottom_friction(kappa,zz,h,z0b,cbcmin,cbcmax);                                                  

%==========================================
%           begin internal (3-D) mode
%==========================================
for iint=1:iend

    [time,ramp] = get_time(iint,time0,lramp,dti,cor);
    
% set time dependent surface boundary conditions     
    [wusurf,wvsurf,wtsurf,wssurf,vfluxf,w] = surface_forcing(iint,time,iend,nsbdy,t,s,vfluxf,wusurf,wvsurf,wtsurf,wssurf,w);

% set time dependent lateral boundary conditions                                                                   
%    lateral_bc                                                                                              
                                                                                                                    
% set lateral viscosity                                                                                            
    [advx,advy,drhox,drhoy,aam] = lateral_viscosity(npg,u,v,dt_3d,aam,aam_init,ub,vb,rho,rmean,ramp,horcon);
                                                                                                                    
% form vertical averages of 3-D fields for use in external (2-D) mode                                              
    [adx2d,ady2d,drx2d,dry2d,aam2d,advua,advva,egf,utf,vtf] = mode_interaction(advx,advy,drhox,drhoy,aam,uab,vab,ua,va,el,d);

%     [e_atmos,vfluxf,wtsurf,swrad,wssurf,advx,advy,drhox,drhoy,aam,adx2d,ady2d,drx2d,    ...
%     dry2d,aam2d,advua,advva,w,egf,utf,vtf]= set_lateralboundary (e_atmos,vfluxf,wtsurf, ...
%     swrad,wssurf,u,ub,v,vb,dt_3d,aam,rho,rmean,ramp,uab,vab,ua,va,d,t,s,w,el);
          
        for iext=1:isplit    % Begin external (2-D) mode

            [elf]=external_el(elb,d,ua,va,vfluxf);

            [uaf,advua] = external_ua(h,elb,el,elf,uab,vab,adx2d,drx2d,advua,d,va,e_atmos,wusurf,wubot,iext,aam2d,ua,ramp);
            
            [vaf,advva] = external_va(h,elb,el,elf,uab,vab,ady2d,dry2d,advva,d,ua,e_atmos,wvsurf,wvbot,iext,aam2d,va,ramp);

            [d,d_3d,ua,va,el,uab,vab,elb,egf,vtf,utf,etf,vamax,imax,jmax]=external_update(uaf,vaf,elf,ua,va,el,uab,...
                                                                                          vab,elb,iext,egf,vtf,utf,etf);
        end

        %===========================================
        %End of external (2-D) mode
        %=============================================        
        
    if(vamax<=vmaxl) 

        if((iint~= 1|| time0~=0.e0) && mode~=2)

            [u,v] = adjust_uv(u,utb,utf,v,vtb,vtf,dt_3d);
          
            w     = internal_vertvl (dt_3d,u,v,vfluxb,vfluxf,etf,etb);  

            [q2f,q2,q2b,q2lf,q2l,q2lb,km,kq,kh]=internal_q(q2b,q2,q2lb,q2l,dt_3d,u,v,w,aam,etb,etf,rho,wusurf,wubot,wvsurf,wvbot,km,kq,kh,t,s);

            [tf,t,tb]=internal_t(tb,t,dt,dt_3d,u,v,aam,w,etb,etf,wtsurf,tsurf,nbct,swrad,kh,nadv,tclim);
            
            [sf,s,sb]=internal_s(sb,s,dt,dt_3d,u,v,aam,w,etb,etf,wssurf,ssurf,nbcs,swrad,kh,nadv,sclim); 
            
            [rho]=dens(s,t,h_3d,fsm_3d);
            
            [uf,wubot] =internal_u(advx,dt_3d,e_atmos,drhox,ub,u,v,w,egf,egb,etf,etb,km,vb,wusurf,cbc);
            
            [vf,wvbot] =internal_v(advy,dt_3d,e_atmos,drhoy,vb,u,v,w,egf,egb,etf,etb,km,ub,wvsurf,cbc);
            
            [u,v,ub,vb] = adjust_ufvf(u,ub,uf,v,vb,vf);

        end
        
        [egb,etb,et,dt,utb,vtb,vfluxb,dt_3d]=internal_update(egf,et,etf,utf,vtf,vfluxf);

    end   
        print_section(iint,vamax,time,iext,tb,sb,dt_3d,et,imax,jmax);
    %-----------------------------------------------------------------------
    %  End of internal (3-D) mode
    %-----------------------------------------------------------------------  
end

