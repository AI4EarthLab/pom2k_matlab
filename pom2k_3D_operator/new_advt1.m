function [ff]=new_advt1(fb,f,dt,u,v,aam,tprni,w,etb,etf,h)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is centred scheme, as originally provide in    *
% *                POM (previously called advt).                       *
% *                                                                    *
% **********************************************************************
load('grid.mat');load('para.mat');
xflux=zeros(im,jm,kb); yflux=zeros(im,jm,kb); zflux=zeros(im,jm,kb); ff=zeros(im,jm,kb);
dt_3d=repmat(dt,1,1,kb);
h_3d=repmat(h,1,1,kb);
% for advective fluxes and diffusive fluxes:
xflux=(AXB1(dt_3d).*AXB1(f).*u-DIVISION(AXB1(aam).*AXB1(h_3d)*tprni.*DXB1(fb).*dum_3d , AXB1(dx_3d))).*AXB1(dy_3d);
yflux=(AYB1(dt_3d).*AYB1(f).*v-DIVISION(AYB1(aam).*AYB1(h_3d)*tprni.*DYB1(fb).*dvm_3d , AYB1(dy_3d))).*AYB1(dx_3d);
% for vertical advection:
zflux=AZB1(f).*w.*art_3d;
xflux(:,:,kb)=0.0; yflux(:,:,kb)=0.0; zflux(:,:,kb)=0.0;

%     Add net horizontal fluxes and then step forward in time:
ff= DIVISION((fb.*(h_3d+repmat(etb,1,1,kb)).* art_3d - dti2 .* (DXF1(xflux)+ DYF1(yflux)-DZF1(zflux)./dz_3d) )  ...
        , ((h_3d+repmat(etf,1,1,kb)).* art_3d)) ;
ff(1,:,:)=0;ff(im,:,:)=0; ff(:,1,:)=0;ff(:,jm,:)=0; ff(:,:,kb)=0;
end