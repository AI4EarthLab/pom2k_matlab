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
         [drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp);   
% % %           aam=horcon .* dx_3d .* dy_3d .*sqrt( (DXF(u)./dx_3d).^2 + (DYF(v)./dy_3d).^2    ...
% % %                   +0.5*( DYB(AYF(AXF(u)))./dy_3d + DXB(AXF(AYF(v)))./dx_3d).^2 );     
        aam=horcon .* dx_3d .* dy_3d .*sqrt( DXF(u).^2 + DYF(v).^2+0.5*( DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))) ).^2 );     
        aam(:,1,:)=aam_init;
        aam(:,jm,:)=aam_init;
        aam(1,:,:)=aam_init;
        aam(im,:,:)=aam_init;
        aam(:,:,kb)=aam_init;
