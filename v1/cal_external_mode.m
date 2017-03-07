 %----------------------------------------------------------------------    
    egf=el*ispi;
    utf=ua .* 2.0 .* AXB(d) .* isp2i;
    vtf=va .* 2.0 .* AYB(d) .* isp2i;   
 for iext=1:isplit    % Begin external (2-D) mode        
        %
        %     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
        %     with pom98.f. See also modifications to subroutine vertvl.
        %    
        % compute SSH in the external model equation 
        elf= elb-dte2.*((DXF(AXB(d).*ua)+DYF(AYB(d).*va))-vfluxf);  
        
        [elf,uaf,vaf,uf,vf,w] = bcond(1,elf,uaf,vaf,uf,vf,w,...
            im,jm,kb,imm1,jmm1,kbm1,...
            fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
            dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        
        if(mod(iext,ispadv)==0)
            [advua,advva] = advave(aam2d,uab,vab,ua,va,d);    
        end   
 
        % compute u in the external model equation 
        uaf=    DIVISION((AXB(h+elb) .* uab -2.0* dte .* (adx2d + advua - AXB(cor .* d .* AYF(va)) ...
                     + grav.* AXB(d).*( (1.0-2.0*alpha) .* DXB(el) + alpha* (DXB(elb)+ DXB(elf)) + DXB(e_atmos) ) ...
                     + drx2d +  (wusurf-wubot) ) ), AXB(h+elf)) ;    

        % compute v in the external model equation     
        vaf=   DIVISION((AYB(h+elb) .* vab -2.0* dte .* (ady2d + advva + AYB(cor .* d .* AXF(ua)) ...
                + grav .* AYB(d).*( (1.0-2.0*alpha) .* DYB(el) + alpha* (DYB(elb)+ DYB(elf)) + DYB(e_atmos) ) ...
                + dry2d + (wvsurf-wvbot) )), AYB(h+elf)) ;  


        [elf,uaf,vaf,uf,vf,w] = bcond(2,elf,uaf,vaf,uf,vf,w,...
            im,jm,kb,imm1,jmm1,kbm1,...
            fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
            dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        
        if(iext==(isplit-2))
            etf=0.25*smoth*elf;
        elseif(iext==(isplit-1))
            etf=etf+0.5*(1.0-0.5*smoth)*elf;
        elseif(iext==isplit)
            etf=(etf+0.5*elf).*fsm;
        end
         
        % Stop if velocity condition violated (generally due to CFL
        % criterion not being satisfied):   
        
        [vamax, vapos]=field_max(abs(vaf(:)));
        [imax, jmax]=ind2sub(size(vaf), vapos);

        
        if(vamax<=vmaxl)
            %
            %     Apply filter to remove time split and reset time sequence:
            %
%             [ua,va,uab,vab]=smoth_update(uaf,vaf,ua,va,uab,vab);
            ua=ua+0.5*smoth*(uab-2.0*ua+uaf);
            va=va+0.5*smoth*(vab-2.0*va+vaf);
            uab=ua;
            ua=uaf;
            vab=va;
            va=vaf;
            el=el+0.5*smoth*(elb-2.0*el+elf);
            elb=el;
            el=elf;
            d=h+el;
            d_3d=repmat(d,1,1,kb);  %add by hx
            
            if(iext~=isplit)
                egf=egf+el*ispi;
                utf=utf+2.0* ua .* AXB(d) * isp2i;
                vtf=vtf+2.0* va .* AYB(d) * isp2i;
            end
        end
    end
    
    %===========================================
    %End of external (2-D) mode
    %=============================================