clear all;

global im jm kb imm1 imm2 jmm1 jmm2 kbm1 kbm2 kl1 kl2 ...
    z zz dz dzz z_3d zz_3d dz_3d dzz_3d ...
    dx dy art aru arv  fsm dum dvm ...
    dx_3d dy_3d art_3d aru_3d arv_3d fsm_3d dum_3d dvm_3d;

global small pi netcdf_file iproblem mode nadv nitera ...
     sw nread dte isplit days prtd1 prtd2 swtch iskp jskp ...
     lramp rhoref tbias sbias grav kappa z0b cbcmin cbcmax...
     horcon tprni umol hmax vmaxl slmax ntp nbct nbcs ispadv...
     smoth alpha dti dte2 dti2 iend iprint iswtch ispi isp2i;

init_constants();
init_variables();

if(iproblem ~= 3)
    [z,zz,dz,dzz,z_3d,zz_3d,dz_3d,dzz_3d]=depth(im, jm, kb, kl1, kl2);
end

global OP
OP = gen_operator(im,jm,kb);


if(iproblem == 1)
    [dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs, ...
        dx_3d,dy_3d,cor_3d, ...
        h_3d,art_3d,aru_3d,arv_3d, ...
        fsm_3d,dum_3d,dvm_3d,dt_3d] = seamount(e_atmos,aam);
   
elseif(iproblem == 2)
   [dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,rot,tatm,satm,vfluxf,tbias,sbias] = box(dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,...
        e_atmos,aam,im,jm,kb,imm1,jmm1,kbm1,slmax,zz,tbias,sbias,grav,rhoref,rot,tatm,satm,vfluxf);
     
elseif(iproblem ==3)
[dx,dy,cor,east_c,north_c,east_e,north_e,east_u,north_u,east_v,...
         north_v,h,art,aru,arv,fsm,dum,dvm,tb,sb,...
         tclim,sclim,ub,vb,elb,etb,dt,aam2d,rmean,rfe,rfw,...
         rfn,rfs,uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,...
         tbn,tbs,sbn,sbs,rot,els,eln,vabs,vabn,ubw,...
         ube,vbs,vbn,ssurf,tsurf,dx_3d,dy_3d,cor_3d,h_3d,art_3d,...
         aru_3d,arv_3d,fsm_3d,dum_3d,dvm_3d,dt_3d,z,zz,dz,dzz,...
         z_3d,zz_3d,dz_3d,dzz_3d] = file2ic(aam);
         
else
    fprintf('Invalid value of iproblem ..... program terminated\n');
    return
end

% create grid and setup the grid size
global gs
gs = create_grid(im, jm, kb, 'C');
gs = init_grid(gs, dx_3d, dy_3d, dz_3d);

init_fields();
%     Inertial period for temporal filter:
period=(2.0*pi)/abs(cor(floor(im/2),floor(jm/2)))/86400.0;

%     Initialise time:
time0=0.0;
time=0.0;
%     Initial conditions:
%     NOTE that lateral thermodynamic boundary conditions are often set
%     equal to the initial conditions and are held constant thereafter.
%    Users can of course create variable boundary conditions.

ua=uab;
va=vab;
el=elb;
et=etb;
etf=et;
d=h + el;
dt=h + et;
w(:,:,1)=vfluxf;

d_3d=repmat(d,1,1,kb);
dt_3d=repmat(dt,1,1,kb);


l           = 0.1*dt_3d;
%q2b(:,:,:)  = small;
q2b(:,:,:)  = small;
q2lb        = l.*q2b;
kh          = l.*sqrt(q2b);
km          = kh;
kq          = kh;
aam(:,:,:)  = aam_init;

q2 =q2b;
q2l=q2lb;
t  =tb;
s  =sb;
u  =ub;
v  =vb;

[rho]=dens(s,t,h_3d,fsm_3d);

                      
[rho,drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp);


% % % for k=1:kbm1
% % %     drx2d=drx2d+drhox(:,:,k)*dz(k);
% % %     dry2d=dry2d+drhoy(:,:,k)*dz(k);
% % % end

drx2d=sum(drhox.*dz_3d, 3);
dry2d=sum(drhoy.*dz_3d, 3);

%     Calculate bottom friction coefficient:
cbc=(kappa./log((1.0+zz(kbm1))*h/z0b)).^2;
cbc=max(cbcmin,cbc);
%     If the following is invoked, then it is probable that the wrong
%     choice of z0b or vertical spacing has been made:
cbc=min(cbcmax,cbc);

%     Calculate external (2-D) CFL time step:
tps=0.5./sqrt(1.0./dx.^2+1.0./dy.^2) ./ sqrt(grav*(h+small)) .* fsm;

d = h+el;
dt = h+et;
time = time0;

d_3d=repmat(d,1,1,kb);
dt_3d=repmat(dt,1,1,kb);

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

    %       write(6,2) mode,iint,time
    %   2   format(' mode,iint,time =',2i5,f9.2)
    %-----------------------------------------------------------------------
    %     Set time dependent, surface and lateral boundary conditions.
    %     The latter will be used in subroutine bcond. Users may
    %     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
    %     wssurf, swrad and vflux.
    %
    %     Introduce simple wind stress. Value is negative for westerly or
    %     southerly winds. The following wind stress has been tapered
    %     along the boundary to suppress numerically induced oscilations
    %     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
    %     To make a healthy surface Ekman layer, it would be well to set
    %     kl1=9.
    %
    for j=2:jmm1
        for i=2:imm1
            if(iproblem~=3) % constant wind read in file2ic
    %           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1));
                wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))  ...
                    *0.25*(dvm(i,j+1)+dvm(i-1,j+1)     ...
                    +dvm(i-1,j)+dvm(i,j));
                % --- no wind ----
                %           wusurf(i,j)=0.e0;
                wvsurf(i,j)=0.e0;
            end
            
            e_atmos(i,j)=0.e0;
            vfluxf(i,j)=0.e0;
            %
            %     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
            %     the sea surface. See calculation of elf(i,j) below and subroutines
            %     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
            %     is no net flow across lateral boundaries, the basin volume will be
            %     constant; if also vflux(i,j).ne.0, then, for example, the average
            %     salinity will change and, unrealistically, so will total salt.
            %
            w(i,j,1)=vfluxf(i,j);
            %
            %     Set wtsurf to the sensible heat, the latent heat (which involves
            %     only the evaporative component of vflux) and the long wave
            %     radiation:
            %
            wtsurf(i,j)=0.0;
            %
            %     Set swrad to the short wave radiation:
            %
            swrad(i,j)=0.0;
            %
            %     To account for change in temperature of flow crossing the sea
            %     surface (generally quite small compared to latent heat effect)
            %
            tatm=t(i,j,1)+tbias;    % an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias);
            %
            %     Set the salinity of water vapor/precipitation which enters/leaves
            %     the atmosphere (or e.g., an ice cover)
            %
            satm=0.0  ;
            wssurf(i,j)=  vfluxf(i,j)*(satm-s(i,j,1)-sbias)  ;
            %
        end
    end
    %
    %-----------------------------------------------------------------------
    %
    %     Set lateral viscosity:
    %
    %     If mode=2 then initial values of aam2d are used. If one wishes
    %     to use Smagorinsky lateral viscosity and diffusion for an
    %     external (2-D) mode calculation, then appropiate code can be
    %     adapted from that below and installed just before the end of the
    %     "if(mode.eq.2)" loop in subroutine advave.
    %
    %     %alculate Smagorinsky lateral viscosity:
    %
    %       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
    %                                     +.5*(du/dy+dv/dx)**2) )
    %

    if(mode~=2)
% % %  [advx,advy]=advct(u,v,dx,dy,dt,aam,ub,vb,aru,arv,im,jm,kb,imm1,jmm1,kbm1);
        [advx,advy]=advct(u,v,dt_3d,aam,ub,vb);        
% % %         [rho,drhox,drhoy] = baropg(rho,drhox,drhoy,...
% % %             im,jm,imm1,jmm1,kb,kbm1,grav,...
% % %             zz,dt,dum,dvm,ramp,rmean,dx,dy);   
        [rho,drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp);   
% % %           aam=horcon .* dx_3d .* dy_3d .*sqrt( (DXF(u)./dx_3d).^2 + (DYF(v)./dy_3d).^2    ...
% % %                   +0.5*( DYB(AYF(AXF(u)))./dy_3d + DXB(AXF(AYF(v)))./dx_3d).^2 );     
        aam=horcon .* dx_3d .* dy_3d .*sqrt( DXF(u).^2 + DYF(v).^2    ...
                  +0.5*( DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))) ).^2 );     
        aam(:,1,:)=aam_init;
        aam(:,jm,:)=aam_init;
        aam(1,:,:)=aam_init;
        aam(im,:,:)=aam_init;
        aam(:,:,kb)=aam_init;

        %     Form vertical averages of 3-D fields for use in external (2-D)
        %     mode:
        
        adx2d=zeros(im,jm);
        ady2d=zeros(im,jm);
        drx2d=zeros(im,jm);
        dry2d=zeros(im,jm);
        aam2d=zeros(im,jm);
       %
%         for k=1:kbm1
%             for j=1:jm
%                 for i=1:im
%                     adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k);
%                     ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k);
%                     drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k);
%                     dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k);
%                     aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k);
%                 end
%             end
%         end       
%         for k=1:kb
%             adx2d=adx2d+advx(:,:,k)*dz(k);
%             ady2d=ady2d+advy(:,:,k)*dz(k);
%             drx2d=drx2d+drhox(:,:,k)*dz(k);
%             dry2d=dry2d+drhoy(:,:,k)*dz(k);
%             aam2d=aam2d+aam(:,:,k)*dz(k);
%         end

        adx2d = sum(advx.*dz_3d, 3);
        ady2d = sum(advy.*dz_3d, 3);
        drx2d = sum(drhox.*dz_3d, 3);
        dry2d = sum(drhoy.*dz_3d, 3);
        aam2d = sum(aam.*dz_3d, 3);
        
        %[tps2,advua,advva,fluxua,fluxva,wubot,wvbot,tps0] = advave(tps,advua,advva,fluxua,fluxva,wubot,wvbot,tps,...
        %    mode,im,jm,imm1,jmm1,aam2d,...
        %    uab,vab,dx,dy,ua,va,cbc,aru,arv,d);
        %[tps,wubot,wvbot,advua,advva] = advave(tps,wubot,wvbot,mode,aam2d,uab,vab,ua,va,cbc,d);
        [advua,advva] = advave(aam2d,uab,vab,ua,va,d);
        adx2d=adx2d-advua;
        ady2d=ady2d-advva;
    end
    
    egf=el*ispi;
    
    utf=ua .* 2.0 .* AXB(d) .* isp2i;
    vtf=va .* 2.0 .* AYB(d) .* isp2i;
    
    %----------------------------------------------------------------------    
    for iext=1:isplit    % Begin external (2-D) mode        
        %
        %     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
        %     with pom98.f. See also modifications to subroutine vertvl.
        %   
        %elf= elb+dte2.*(-(DXF( AXB(d).*AXB(dy).*ua)+DYF(AYB(d).*AYB(dx).*va))./ art-vfluxf);  
        
        % compute SSH in the external model equation 
        elf= elb-dte2.*((DXF(AXB(d).*ua)+DYF(AYB(d).*va))-vfluxf);  
        
        [elf,uaf,vaf,uf,vf,w] = bcond(1,elf,uaf,vaf,uf,vf,w,...
            im,jm,kb,imm1,jmm1,kbm1,...
            fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
            dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        
        if(mod(iext,ispadv)==0)
% % %             [tps2,advua,advva,fluxua,fluxva,wubot,wvbot,tps0] = advave(tps,advua,advva,fluxua,fluxva,wubot,wvbot,tps,...
% % %                 mode,im,jm,imm1,jmm1,aam2d,...
% % %                 uab,vab,dx,dy,ua,va,cbc,aru,arv,d);
%%%            [tps,wubot,wvbot,advua,advva] = advave(tps,wubot,wvbot,mode,aam2d,uab,vab,ua,va,cbc,d); 
            [advua,advva] = advave(aam2d,uab,vab,ua,va,d);    
        end
    
% % %         uaf=   DIVISION( (2.0*AXB(h+elb) .* aru .* uab ...
% % %                -4.0* dte .* (adx2d + advua - aru .* AXB(cor .* d .* AYF(va)) ...
% % %                + grav .* AXB(dy) .* AXB(d)  ...
% % %                  .*( (1.0-2.0*alpha) .* DXB(el) + alpha* (DXB(elb)+ DXB(elf)) + DXB(e_atmos) ) ...
% % %                + drx2d + aru .* (wusurf-wubot))) , (2.0*AXB(h+elf) .* aru));      
 
        % compute u in the external model equation 
        uaf=    (AXB(h+elb) .* uab ...
                    -2.0* dte .*  ...
                    (adx2d + advua - AXB(cor .* d .* AYF(va)) ...
                     + grav.* AXB(d).*( (1.0-2.0*alpha) .* DXB(el) + alpha* (DXB(elb)+ DXB(elf)) + DXB(e_atmos) ) ...
                     + drx2d +  (wusurf-wubot) ...
                    ) ...
                 )./ AXB(h+elf) ;    

% % %         vaf=   DIVISION( (2.0*AYB(h+elb) .* arv .* vab ...
% % %                -4.0* dte .* (ady2d + advva + arv .* AYB(cor .* d .* AXF(ua)) ...
% % %                + grav .* AYB(dx) .* AYB(d)  ...
% % %                  .*( (1.0-2.0*alpha) .* DYB(el) + alpha* (DYB(elb)+ DYB(elf)) + DYB(e_atmos) ) ...
% % %                + dry2d + arv .* (wvsurf-wvbot))) , (2.0*AYB(h+elf) .* arv));  

        % compute v in the external model equation     
        vaf=   (AYB(h+elb) .* vab ...
                    -2.0* dte .* ...
                    (ady2d + advva + AYB(cor .* d .* AXF(ua)) ...
                     + grav .* AYB(d).*( (1.0-2.0*alpha) .* DYB(el) + alpha* (DYB(elb)+ DYB(elf)) + DYB(e_atmos) ) ...
                     + dry2d + (wvsurf-wvbot) ...
                    ) ...
               )./ AYB(h+elf) ;  


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
        
%        vamax=0.0;       
%         for j=1:jm
%             for i=1:im
%                 if(abs(vaf(i,j))>=vamax)
%                     vamax=abs(vaf(i,j));
%                     imax=i;
%                     jmax=j;
%                 end
%             end
%         end
        
        [vamax, vapos]=field_max(abs(vaf(:)));
        [imax, jmax]=ind2sub(size(vaf), vapos);
%         fprintf('diff_vmax:%.18f\n', vamax-vamax);
%         fprintf('diff_imax:%d\n', imax1-imax);
%         fprintf('diff_jmax:%d\n', jmax1-jmax);
        
        if(vamax<=vmaxl)
            %
            %     Apply filter to remove time split and reset time sequence:
            %
            ua=ua+0.5*smoth*(uab-2.0*ua+uaf);
            va=va+0.5*smoth*(vab-2.0*va+vaf);
            el=el+0.5*smoth*(elb-2.0*el+elf);
            elb=el;
            el=elf;
            d=h+el;
            uab=ua;
            ua=uaf;
            vab=va;
            va=vaf;
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
    
    if(vamax<=vmaxl) 
        %
        %     continue with internal (3-D) mode calculation:
        %
        if((iint~= 1|| time0~=0.e0) && mode~=2)
            %
            %     Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
            %
%             tps=zeros(im,jm);
%             for k=1:kbm1
%                 tps=tps+u(:,:,k)*dz(k);
%             end
%?
            tps=sum(u(:,:,1:kbm1).*dz_3d(:,:,1:kbm1), 3);
            tps1=sum(u.*dz_3d, 3);            
%             for k=1:kbm1
%                 u(:,:,k)=(u(:,:,k)-tps)+   DIVISION((utb+utf), 2.0*AXB(dt));
%             end
            
            utb_3d = repmat(utb, 1, 1, kbm1);
            utf_3d = repmat(utf, 1, 1, kbm1);
            tps_3d = repmat(tps, 1, 1, kbm1);
            dt_3d1 = repmat(dt, 1, 1, kbm1);
            dt_axb = 2.0 * AXB(dt_3d1);
            u(:,:,1:kbm1) = u(:,:,1:kbm1)-tps_3d + (utb_3d+utf_3d) ./ dt_axb;
            
% 
%             tps =zeros(im,jm);
%             for k=1:kbm1
%                 tps=tps+v(:,:,k)*dz(k);
%             end
%?
            tps = sum(v(:,:,1:kbm1) .* dz_3d(:,:,1:kbm1), 3);
            
%             for k=1:kbm1
%                 v(:,:,k)=(v(:,:,k)-tps) + DIVISION((vtb+vtf), 2.0*AYB(dt));
%             end
           
            vtb_3d = repmat(vtb, 1, 1, kbm1);
            vtf_3d = repmat(vtf, 1, 1, kbm1);
            tps_3d = repmat(tps, 1, 1, kbm1);
            dt_ayb = 2.0 * AYB(dt_3d1);
            v(:,:,1:kbm1) = v(:,:,1:kbm1) - tps_3d + (vtb_3d+vtf_3d) ./ dt_ayb;
            

            %     vertvl calculates w from u, v, dt (h+et), etf and etb:
            %
            dt_3d=repmat(dt,1,1,kb);
            [a,c,w]=vertvl (w,dt_3d,u,v,vfluxb,vfluxf,etf,etb,dti2);  
            %[a0,c0,w0]=vertvl(a,c,w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1);  

            [elf,uaf,vaf,uf,vf,w] = bcond(5,elf,uaf,vaf,uf,vf,w,...
                im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
            
            vf=zeros(im,jm,kb);
            uf=zeros(im,jm,kb);
                          
            % calculate q2f and q2lf using uf, vf, a and c as temporary variables:
           % [q2b,q2,uf,a,c]=advq(q2b,q2,uf,a,c,...
           %     dt,dx,dy,dz,u,v,w,aam,h,dum,dvm,art,etb,etf,im,jm,imm1,jmm1,kbm1,dti2);
           uf=advq(q2b,q2,dt_3d,u,v,w,aam,h,etb,etf,dti2);    
           % [q2lb,q2l,vf,a,c]=advq(q2lb,q2l,vf,a,c,...
           %     dt,dx,dy,dz,u,v,w,aam,h,dum,dvm,art,etb,etf,im,jm,imm1,jmm1,kbm1,dti2);
           vf=advq(q2lb,q2l,dt_3d,u,v,w,aam,h,etb,etf,dti2);  
          
%             [a,c,tps,dtef,...
%                 ee,gg,l,kq,km,kh,...
%                 uf,vf,q2b,q2lb,a,c]=profq(a,c,tps,dtef,....
%                 ee,gg,l,kq,km,kh,...
%                 uf,vf,q2,q2b,q2lb,a,c,...
%                 h,etf,dti2,umol,dzz,grav,rho,kappa,u,v,dt,small,fsm,im,jm,kb,imm1,jmm1,kbm1,tbias,sbias,dz,...
%                 wusurf,wubot,wvsurf,wvbot,t,s,rhoref,zz,z);

            [sm,sh,dh,l,kq,km,kh,uf,vf,q2b,q2lb]=profq(l,kq,km,kh,uf,...
           vf,q2,q2b,q2lb,etf,h,rho,u,v,dt,wusurf,wubot,wvsurf,wvbot,t,s);
            
            [elf,uaf,vaf,uf,vf,w] = bcond(6,elf,uaf,vaf,uf,vf,w,...
                im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
            
            q2=q2+.5e0*smoth*(uf+q2b-2.e0*q2);
            q2l=q2l+.5e0*smoth*(vf+q2lb-2.e0*q2l);
            q2b=q2;
            q2=uf;
            q2lb=q2l;
            q2l=vf;
            
            %
            %      calculate tf and sf using uf, vf, a and c as temporary variables:
            %
            if(mode~=4)
                if(nadv==1)    
%                     uf=advt1(tb,t,dt,u,v,aam,tprni,w,etb,etf,h);
                    uf=advt1(tb,t,dt_3d,u,v,aam,w,etb,etf,h_3d);
%                     vf=advt1(sb,s,dt,u,v,aam,tprni,w,etb,etf,h);
                    vf=advt1(sb,s,dt_3d,u,v,aam,w,etb,etf,h_3d);
                elseif(nadv==2)                  
                    [tb,t,tclim,uf,a,c,nitera,sw,...
                     zflux] = advt2(tb,t,tclim,uf,a,c,nitera,sw,...
                                                 zflux,...
                                                 im,jm,kb,imm1,jmm1,kbm1,dti2,...
                                                 etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);
                    [sb,s,sclim,vf,a,c,nitera,sw,...
                     zflux] = advt2(sb,s,sclim,vf,a,c,nitera,sw,...
                                                 zflux,...
                                                 im,jm,kb,imm1,jmm1,kbm1,dti2,...
                                                 etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);      
                    %
                else
                    disp 'Invalid value for nadv ..... '
                    disp 'program terminated'
                    
                    return
                end

                [uf] = proft(uf,wtsurf,tsurf,nbct,h,etf,swrad,kh);
                [vf] = proft(vf,wssurf,ssurf,nbcs,h,etf,swrad,kh);
                
                [elf,uaf,vaf,uf,vf,w] = bcond(4,elf,uaf,vaf,uf,vf,w,...
                    im,jm,kb,imm1,jmm1,kbm1,...
                    fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                    dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
                               
                t=t+.5e0*smoth*(uf+tb-2.e0*t);
                s=s+.5e0*smoth*(vf+sb-2.e0*s);
                tb=t;
                t=uf;
                sb=s;
                s=vf;

                [rho]=dens(s,t,h_3d,fsm_3d);
            end  % end if
            % calculate uf and vf:
%           uf = advu(advx,cor,dt_3d,e_atmos,drhox,h_3d,ub,u,v,w,egf,egb,etf,etb);
            uf = advu(advx,cor,dt_3d,e_atmos,drhox,h_3d,ub,u,v,w,egf,egb,etf,etb);
            vf = advv(advy,cor,dt_3d,e_atmos,drhoy,h_3d,vb,u,v,w,egf,egb,etf,etb);
            
            vf0=vf;
            
            [uf,wubot] = profu(uf,etf,h,km,wusurf,cbc,ub,vb);           
            [vf,wvbot] = profv(vf,etf,h,km,wvsurf,cbc,ub,vb);
            
            [elf,uaf,vaf,uf,vf,w] = bcond(3,elf,uaf,vaf,uf,vf,w,...
                im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
            
%             tps = zeros(im,jm);
% 
%             for k=1:kbm1
%                 for j=1:jm
%                     for i=1:im
%                         tps(i,j)=tps(i,j)     ...
%                             +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k);
%                     end
%                 end
%             end
            tps = sum((uf+ub-2.e0*u).*dz_3d, 3);
            %
%             for k=1:kbm1
%                 for j=1:jm
%                     for i=1:im
%                         u(i,j,k)=u(i,j,k)     ...
%                             +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)     ...
%                             -2.e0*u(i,j,k)-tps(i,j));
%                     end
%                 end
%             end
            tps_3d = repmat(tps, 1, 1, kbm1);
            u(:,:,1:kbm1)=u(:,:,1:kbm1)+ ...
                .5e0*smoth*(uf(:,:,1:kbm1)+ub(:,:,1:kbm1)-2.e0*u(:,:,1:kbm1)-tps_3d);
            %
%           tps = zeros(im,jm);
%             for k=1:kbm1
%                 for j=1:jm
%                     for i=1:im
%                         tps(i,j)=tps(i,j)     ...
%                             +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k);
%                     end
%                 end
%             end
            tps = sum((vf+vb-2.e0*v).*dz_3d, 3);
            
            %
%             for k=1:kbm1
%                 for j=1:jm
%                     for i=1:im
%                         v(i,j,k)=v(i,j,k)     ...
%                             +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)     ...
%                             -2.e0*v(i,j,k)-tps(i,j));
%                     end
%                 end
%             end
            %
            tps_3d = repmat(tps,1,1,kbm1);
            v(:,:,1:kbm1) = v(:,:,1:kbm1)+ ...
                .5e0*smoth*(vf(:,:,1:kbm1)+vb(:,:,1:kbm1)-2.e0*v(:,:,1:kbm1)-tps_3d);
            
            ub=u;
            u=uf;
            vb=v;
            v=vf;
            %
        end  % end if
        %        %
        egb=egf;
        etb=et;
        et=etf;
        dt=h+et;
        utb=utf;
        vtb=vtf;
        vfluxb=vfluxf;
        dt_3d=repmat(dt,1,1,kb);        %add by hx
	%
    end   % end if
    
  

    %
    %     Beginning of print section:
    %
    if(iint>=iswtch)
        iprint=floor(prtd2*24.e0*3600.e0/dti + 0.5);
    end
    %
    if(mod(iint,iprint)==0 || vamax>=vmaxl)
        fprintf('**************************************************\n');
        fprintf('**************************************************\n');
        fprintf('time = %.6f, iint = %d, iexit = %d, iprint = %d\n',time,iint,iext,iprint);
        
        %
        %     Select print statements in printall as desired:
        %
        %printall
        %
        vtot=0.e0;
        atot=0.e0;
        taver=0.e0;
        saver=0.e0;
        eaver=0.e0;
%         for k=1:kbm1
%             for j=1:jm
%                 for i=1:im
%                     darea=dx(i,j)*dy(i,j)*fsm(i,j);
%                     dvol=darea*dt(i,j)*dz(k);
%                     vtot=vtot+dvol;
%                     taver=taver+tb(i,j,k)*dvol;
%                     saver=saver+sb(i,j,k)*dvol;
%                 end
%             end
%         end
        
        dt_3d1 = repmat(dt, 1, 1, kb);
        tmp_dvol = dx_3d(:,:,1:kbm1).*dy_3d(:,:,1:kbm1).*fsm_3d(:,:,1:kbm1).*dt_3d1(:,:,1:kbm1).*dz_3d(:,:,1:kbm1);
        vtot=sum(tmp_dvol(:));
        taver=sum(reshape(tb(:,:,1:kbm1).*tmp_dvol,[],1));
        saver=sum(reshape(sb(:,:,1:kbm1).*tmp_dvol,[],1));

%         fprintf('diff_vtot:%.18f\n', vtot1-vtot);
%         fprintf('diff_taver:%.18f\n', taver1-taver);
%         fprintf('diff_saver:%.18f\n', saver1-saver);
        
        %
%         for j=1:jm
%             for i=1:im
%                 darea=dx(i,j)*dy(i,j)*fsm(i,j);
%                 atot=atot+darea;
%                 eaver=eaver+et(i,j)*darea;
%             end
%         end
        %
        darea = dx.*dy.*fsm;
        atot = sum(darea(:));
        eaver = sum(reshape(et.*darea,[],1));
        
        taver=taver/vtot;
        saver=saver/vtot;
        eaver=eaver/atot;
        tsalt=(saver+sbias)*vtot;
      
        fprintf('vtot = %.6f,atot = %.6f\n',vtot,atot);
        fprintf('eaver = %.6f,taver = %.6f,saver=%.6f,saver = %.6f,tsalt = %.6f\n',eaver,taver,saver,tsalt,tsalt);
        
        
        if(vamax>vmaxl)
            
            fprintf('time = %.6f, iint = %.6f, iexit = %.6f, iprint = %.6f\n',time,iint,iext,iprint);
           
%             printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s,rho,aam,km,kb,mode,dt,zz,z);
           
            disp        '************************************************'
            disp        '************ abnormal job end ******************'
            disp        '************* user terminated ******************'
            disp        '************************************************'
            
            fprintf('vamax = %d, imax = %d ,jmax = %d \n',vamax,imax,jmax);
            return;
            %
        end
        %
    end
    %
    %     End of print section
    %-----------------------------------------------------------------------
    %
    %  End of internal (3-D) mode
    %-----------------------------------------------------------------------  
end
%internal loop ended
%write(6,4) time,iint,iext,iprint

fprintf('time = %d, iint = %d ,iext = %d , iprint = %d \n',time,iint,iext,iprint);
%
%     Set levels for output:
%
ko(1)=1;
ko(2)=2;
ko(3)=kb/2;
ko(4)=kb-1;
ko(5)=kb;
%
%      prxyz('Vertical velocity, w                    ',
%    ...           time,w       ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
%
%      prxyz('Turbulent kinetic energy x 2, q2        ',
%    ...           time,q2      ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
%
%     Save this data for a seamless restart:
%
%    write_to_file(im,jm,kb,time,wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,...
%    utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,...
%    km,kh,kq,l,q2,q2b,aam,q2l,q2lb);
%    printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s,rho,aam,km,kb,mode,dt,zz,z);
    
%     End of main program
