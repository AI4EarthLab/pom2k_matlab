clear all;
init_global();
init_constants();
init_operator();
create_grid();
init_variables();
read_initvalues();
init_grid();
init_fields();

[period,time,time0,ua,va,el,et,etf,d,dt,w,d_3d,dt_3d,l,q2b,q2lb,kh,km,kq,aam,   ...
 q2,q2l,t,s,u,v,rho,drhox,drhoy,drx2d,dry2d,cbc,tps] = init_time_condition(uab, ...
 vab,elb,etb,vfluxf,tb,sb,ub,vb,rmean,ramp,q2b,aam,w);

%==========================================
%           begin internal (3-D) mode
%==========================================
for iint=1:iend
    time=dti*iint*1.0/86400.0+time0;
    if(lramp~=0)
        ramp = time/period;
        if(ramp>1.0)
            ramp=1.0;
        end
    else
        ramp=1.0;
    end
  
    [e_atmos,vfluxf,wtsurf,swrad,wssurf,advx,advy,drhox,drhoy,aam,adx2d,ady2d,drx2d,    ...
    dry2d,aam2d,advua,advva,w,egf,utf,vtf]= set_lateralboundary (e_atmos,vfluxf,wtsurf, ...
    swrad,wssurf,u,ub,v,vb,dt_3d,aam,rho,rmean,ramp,uab,vab,ua,va,d,t,s,w,el);
          
        for iext=1:isplit    % Begin external (2-D) mode

            [elf]=cal_el(elb,d,ua,va,vfluxf);

            [uaf,advua] = cal_ua(h,elb,el,elf,uab,vab,adx2d,drx2d,advua,d,va,e_atmos,wusurf,wubot,iext,aam2d,ua);
            [vaf,advva] = cal_va(h,elb,el,elf,uab,vab,ady2d,dry2d,advva,d,ua,e_atmos,wvsurf,wvbot,iext,aam2d,va);
            [uaf, vaf] = bcond2(uaf,uabe, uabw, rfe, rfw, vaf, vabn,vabs, rfn, rfs, el, ele, eln, elw, els,ramp);

            [d,d_3d,ua,va,el,uab,vab,elb,egf,vtf,utf,etf,vamax,imax,jmax]=external_update(uaf,vaf,elf,ua,...
                    va,el,uab,vab,elb,iext,egf,vtf,utf,etf);
        end

        %===========================================
        %End of external (2-D) mode
        %=============================================        
        
    if(vamax<=vmaxl) 

        if((iint~= 1|| time0~=0.e0) && mode~=2)

            [u,v] = adjust_uv(u,utb,utf,v,vtb,vtf,dt_3d);
          
            w     = vertvl (dt_3d,u,v,vfluxb,vfluxf,etf,etb);  

            [q2f,q2,q2b,q2lf,q2l,q2lb,km,kq,kh]=cal_q(q2b,q2,q2lb,q2l,dt_3d,u,...
                    v,w,aam,etb,etf,rho,wusurf,wubot,wvsurf,wvbot,km,kq,kh,t,s);

            [tf]=cal_ts(tb,t,dt_3d,u,v,aam,w,etb,etf,wtsurf,tsurf,nbct,swrad,kh);
            [sf]=cal_ts(sb,s,dt_3d,u,v,aam,w,etb,etf,wssurf,ssurf,nbcs,swrad,kh); 
            [tf,sf,t,s,tb,sb] = bcond4(tf,sf,tb,sb,u,v,w,t,s,dt,tbe,tbw,tbn,tbs,sbe,sbw,sbn,sbs);
            
            [rho]=dens(s,t,h_3d,fsm_3d);
            
            [uf,wubot] =cal_u(advx,dt_3d,e_atmos,drhox,ub,u,v,w,egf,egb,etf,etb,km,vb,wusurf,cbc);
            [vf,wvbot] =cal_v(advy,dt_3d,e_atmos,drhoy,vb,u,v,w,egf,egb,etf,etb,km,ub,wvsurf,cbc);
            [uf, vf] = bcond3(uf, u, vf, v);
            [u,v,ub,vb] = adjust_ufvf(u,ub,uf,v,vb,vf);

        end
        
        [egb,etb,et,dt,utb,vtb,vfluxb,dt_3d]=internal_update(egf,et,etf,utf,vtf,vfluxf);

    end   
        print_section(iint,vamax,time,iext,tb,sb,dt_3d,et,imax,jmax);
    %-----------------------------------------------------------------------
    %  End of internal (3-D) mode
    %-----------------------------------------------------------------------  
end

