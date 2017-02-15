function [xflux,yflux,w]=vertvl (w,dt_3d,u,v,vfluxb,vfluxf,etf,etb,dti2)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
%load('grid.mat'); load('operator.mat');
global  im kb dz_3d;
%     Reestablish boundary conditions:
% % % xflux = AXB(dy_3d) .* AXB(dt_3d).* u;
% % % yflux = AYB(dx_3d) .* AYB(dt_3d).* v;
xflux = AXB(dt_3d).* u;
yflux = AYB(dt_3d).* v;
%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
w(:,:,1)=0.5e0*(vfluxb+vfluxf);
w=repmat(w(:,:,1),1,1,kb);              %by hx
 
tps=(etf-etb)/dti2;
tps=repmat(tps,1,1,kb);
w=w+SUM2(dz_3d .* ( ( DXF(AXB(dt_3d).* u)+DYF(yflux) )+ tps));  %by hx
w(1,:,:) = 0.0; w(im,:,:) = 0.0;

end