function [xflux,yflux,w]=vertvl (xflux,yflux,w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
%load('grid.mat'); load('operator.mat');
global  kb dx_3d dy_3d dz_3d;
dt_3d=repmat(dt,1,1,kb);
%     Reestablish boundary conditions:
xflux = AXB(dy_3d) .* AXB(dt_3d).* u;
yflux = AYB(dx_3d) .* AYB(dt_3d).* v;
%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
w(:,:,1)=0.5e0*(vfluxb+vfluxf);
w=repmat(w(:,:,1),1,1,kb);              %by hx
 
tps=(etf-etb)/dti2;
tps=repmat(tps,1,1,kb);
w=w+SUM2(dz_3d .* ( ( DXF(xflux)*1.e0+DYF(yflux) )./(dx_3d .*dy_3d)+ tps));  %by hx
w(1,:,:) = 0.e0; w(im,:,:) = 0.e0;
return