function [xflux,yflux,w]=new_vertvl (xflux,yflux,w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('operator.mat');
dt_3d=repmat(dt,1,1,kb);
w = zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);

%     Reestablish boundary conditions:
xflux = AXB(dy_3d) .* AXB(dt_3d) .* u; 
yflux = AYB(dx_3d) .* AYB(dt_3d) .* v;

%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
%
w(:,:,1)=OP_L_XY * (0.5*(vfluxb+vfluxf)) *OP_R_XY;
w0=repmat(w(:,:,1),1,1,kb);

tps=OP_L_XZ*(etf-etb)/dti2;
tps=repmat(tps,1,1,kb);
w = SUM2(dz_3d .* ( ( DXF(xflux)+DYF(yflux) )./(dx_3d.*dy_3d)+ tps))+w0;
w(1,:,:) = 0.e0;

return
