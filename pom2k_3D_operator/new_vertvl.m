function [xflux,yflux,w]=new_vertvl (w,dt,u,v,vfluxb,vfluxf,etf,etb,dti2)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  calculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('operator.mat');

dt_3d=repmat(dt,1,1,kb);
w = zeros(im,jm,kb);
w(:,:,1)=0.5*(vfluxb+vfluxf);
%     Reestablish boundary conditions:      
xflux = AXB1(dy_3d) .* AXB1(dt_3d) .* u; 
yflux = AYB1(dx_3d) .* AYB1(dt_3d) .* v;
% NOTE that, if one wishes to include freshwater flux, the
% surface velocity should be set to vflux(i,j). See also
% change made to 2-D volume conservation equation which
% calculates elf.
tps=OP_L_XZ*(etf-etb)/dti2;
tps=repmat(tps,1,1,kb);
w=SUM2(w+dz_3d .* ( ( DXF2(xflux)+DYF2(yflux) )./(dx_3d.*dy_3d)+ tps));
end