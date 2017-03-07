            uf = advu(advx,cor,dt_3d,e_atmos,drhox,ub,u,v,w,egf,egb,etf,etb);
            vf = advv(advy,cor,dt_3d,e_atmos,drhoy,vb,u,v,w,egf,egb,etf,etb);
                        
            [uf,wubot] = profu(uf,etf,h,km,wusurf,cbc,ub,vb);           
            [vf,wvbot] = profv(vf,etf,h,km,wvsurf,cbc,ub,vb);
            
            [elf,uaf,vaf,uf,vf,w] = bcond(3,elf,uaf,vaf,uf,vf,w,im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);

            adjust_ufvf();