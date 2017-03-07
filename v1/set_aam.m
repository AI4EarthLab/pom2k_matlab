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
        aam=horcon .* dx_3d .* dy_3d .*sqrt( DXF(u).^2 + DYF(v).^2 +0.5*( DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))) ).^2 );     
        aam(:,1,:)=aam_init;
        aam(:,jm,:)=aam_init;
        aam(1,:,:)=aam_init;
        aam(im,:,:)=aam_init;
        aam(:,:,kb)=aam_init;