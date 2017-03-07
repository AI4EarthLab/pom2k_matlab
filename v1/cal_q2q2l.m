           q2f     = advq (q2b,q2,dt_3d,u,v,w,aam,etb,etf);    
           q2lf    = advq (q2lb,q2l,dt_3d,u,v,w,aam,etb,etf);  

            [l,kq,km,kh,q2f,q2lf,q2b,q2lb]=profq(kq,km,kh,q2f,...
           q2lf,q2,q2b,q2lb,etf,rho,u,v,dt,wusurf,wubot,wvsurf,wvbot,t,s);
            
            [elf,uaf,vaf,q2f,q2lf,w] = bcond(6,elf,uaf,vaf,q2f,q2lf,w,im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);

            [q2,q2l,q2b,q2lb]=smoth_update(q2f,q2lf,q2,q2l,q2b,q2lb); 