function [xflux,yflux,w]=new_vertvl (w,dt,u,v,vfluxb,vfluxf,etf,etb,dti2)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  %alculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('operator.mat');
xflux = zeros(im,jm,kb);
yflux = zeros(im,jm,kb);
w = zeros(im,jm,kb);
w(:,:,1)=0.5*(vfluxb+vfluxf);
%     Reestablish boundary conditions:  
for k=1:kbm1
    xflux(:,:,k) = AXB1_XY(dy) .* AXB1_XY(dt) .* u(:,:,k); 
    yflux(:,:,k) = AYB1_XY(dx) .* AYB1_XY(dt) .* v(:,:,k);
    % NOTE that, if one wishes to include freshwater flux, the
    % surface velocity should be set to vflux(i,j). See also
    % change made to 2-D volume conservation equation which
    % calculates elf.
    w(:,:,k+1)=w(:,:,k)+repmat(dz(k),im,jm) .* ( ( DXF2_XY(xflux(:,:,k))+DYF2_XY(yflux(:,:,k)) )./(dx.*dy)+OP_L_XZ*(etf-etb)/dti2 );
end
