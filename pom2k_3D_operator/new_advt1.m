function [ff]=new_advt1(fb,f,dt,u,v,aam,tprni,w,etb,etf,h)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is centred scheme, as originally provide in    *
% *                POM (previously called advt).                       *
% *                                                                    *
% **********************************************************************
%load('grid.mat');load('para.mat');
global im jm kb dum_3d dvm_3d dx_3d dy_3d dz_3d art_3d dti2;

xflux=zeros(im,jm,kb); yflux=zeros(im,jm,kb); zflux=zeros(im,jm,kb); ff=zeros(im,jm,kb);
dt_3d=repmat(dt,1,1,kb);
h_3d=repmat(h,1,1,kb);
% for advective fluxes and diffusive fluxes:
xflux=(AXB(dt_3d).*AXB(f).*u-DIVISION(AXB(aam).*AXB(h_3d)*tprni.*DXB(fb).*dum_3d , AXB(dx_3d))).*AXB(dy_3d);
yflux=(AYB(dt_3d).*AYB(f).*v-DIVISION(AYB(aam).*AYB(h_3d)*tprni.*DYB(fb).*dvm_3d , AYB(dy_3d))).*AYB(dx_3d);
% for vertical advection:
zflux=AZB(f).*w.*art_3d;
xflux(:,:,kb)=0.0; yflux(:,:,kb)=0.0; zflux(:,:,kb)=0.0;

%     Add net horizontal fluxes and then step forward in time:
ff= DIVISION((fb.*(h_3d+repmat(etb,1,1,kb)).* art_3d - dti2 .* (DXF(xflux)+ DYF(yflux)-DZF(zflux)./dz_3d) )  ...
        , ((h_3d+repmat(etf,1,1,kb)).* art_3d)) ;
ff(1,:,:)=0;ff(im,:,:)=0; ff(:,1,:)=0;ff(:,jm,:)=0; ff(:,:,kb)=0;
end