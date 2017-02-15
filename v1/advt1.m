function [ff]=advt1(fb,f,dt_3d,u,v,aam,w,etb,etf,h_3d)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is centred scheme, as originally provide in    *
% *                POM (previously called advt).                       *
% *                                                                    *
% **********************************************************************
%load('grid.mat');load('para.mat');
global im jm kb dum_3d dvm_3d gs dti2 tprni;

xflux=create_field(zeros(im,jm,kb),gs,2); 
yflux=create_field(zeros(im,jm,kb),gs,1); 
zflux=create_field(zeros(im,jm,kb),gs,7); 
ff=create_field(zeros(im,jm,kb),gs,3); 

% for advective fluxes and diffusive fluxes:
% xflux=(AXB(dt_3d).*AXB(f).*u-DIVISION(AXB(aam).*AXB(h_3d)*tprni.*DXB(fb).*dum_3d , AXB(dx_3d))).*AXB(dy_3d);
xflux=AXB(dt_3d).*AXB(f).*u-AXB(aam).*AXB(h_3d)*tprni.*DXB(fb).*dum_3d;
% yflux=(AYB(dt_3d).*AYB(f).*v-DIVISION(AYB(aam).*AYB(h_3d)*tprni.*DYB(fb).*dvm_3d , AYB(dy_3d))).*AYB(dx_3d);
yflux=AYB(dt_3d).*AYB(f).*v-AYB(aam).*AYB(h_3d)*tprni.*DYB(fb).*dvm_3d;
% for vertical advection:
zflux=AZB(f).*w;
xflux(:,:,kb)=0.0; yflux(:,:,kb)=0.0; zflux(:,:,kb)=0.0;

%     Add net horizontal fluxes and then step forward in time:
% ff= DIVISION((fb.*(h_3d+repmat(etb,1,1,kb)).* art_3d - dti2 .* (DXF(xflux)+ DYF(yflux)-DZF(zflux)./dz_3d) )  ...
%         , ((h_3d+repmat(etf,1,1,kb)).* art_3d)) ;
ff= ((h_3d+repmat(etb,1,1,kb)).*fb - dti2 .* (DXF(xflux)+ DYF(yflux)-DZF(zflux)) )  ...
        ./((h_3d+repmat(etf,1,1,kb))) ;        
ff(1,:,:)=0;ff(im,:,:)=0; ff(:,1,:)=0;ff(:,jm,:)=0; ff(:,:,kb)=0;
end