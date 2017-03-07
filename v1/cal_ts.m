            if(mode~=4)
                if(nadv==1)    
                    tf=advt1(tb,t,dt_3d,u,v,aam,w,etb,etf);
                    sf=advt1(sb,s,dt_3d,u,v,aam,w,etb,etf);
                elseif(nadv==2)                  
                    [tb,t,tclim,tf,a,c,nitera,sw,zflux] = advt2(tb,t,tclim,tf,a,c,nitera,sw,zflux,im,jm,kb,imm1,jmm1,kbm1,dti2,...
                                    etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);
                    [sb,s,sclim,sf,a,c,nitera,sw,zflux] = advt2(sb,s,sclim,sf,a,c,nitera,sw,zflux,im,jm,kb,...
                                                 imm1,jmm1,kbm1,dti2,etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);      

                else
                    disp 'Invalid value for nadv ..... '
                    disp 'program terminated'
                    
                    return
                end
                
                [tf] = proft(tf,wtsurf,tsurf,nbct,etf,swrad,kh);
                [sf] = proft(sf,wssurf,ssurf,nbcs,etf,swrad,kh);
                
                [elf,uaf,vaf,tf,sf,w] = bcond(4,elf,uaf,vaf,tf,sf,w,im,jm,kb,imm1,jmm1,kbm1,...
                    fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                    dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);

                [t,s,tb,sb]=smoth_update(tf,sf,t,s,tb,sb);

                [rho]=dens(s,t,h_3d,fsm_3d);
            end  % end if