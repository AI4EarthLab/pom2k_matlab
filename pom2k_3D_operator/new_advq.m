function [qf]=new_advq(qb,q,dt,u,v,w,aam,h,etb,etf,dti2)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates horizontal advection and diffusion, and  *
% *                vertical advection for turbulent quantities.        *
% *                                                                    *
% **********************************************************************
load('grid.mat');load('operator.mat');
xflux = zeros(im,jm,kb);
yflux = zeros(im,jm,kb);
qf=zeros(im,jm,kb);

dt_3d=repmat(dt,1,1,kb);
h_3d =repmat(h,1,1,kb);

%     Do horizontal advection and diffusion:
xflux=AXB1(dy_3d) .* ( AXB1(q) .* AXB1(dt_3d) .* AZB2(u) -AZB2( AXB1( aam ))...
         .*AXB1(h_3d) .*DXB1( qb ) .* dum_3d .* DIVISION( 1.0,AXB1( dx_3d ) ) );
%     do vertical advection: add flux terms, then step forward in time:
yflux=AYB1(dx_3d) .* ( AYB1(q) .* AYB1(dt_3d) .* AZB2(v) -AZB2( AYB1( aam ) )...
         .*AYB1(h_3d) .*DYB1( qb ) .* dvm_3d .* DIVISION( 1.0,AYB1( dy_3d ) ) );

qf= ( ( h_3d+repmat(etb,1,1,kb) ) .* art_3d .* qb...
          -dti2*( -DZC(w.*q) .* art_3d.*DIVISION( 1.0,AZB2(dz_3d) ) + DXF2(xflux)+DYF2(yflux)) )...
          ./ ( ( h_3d+repmat(etf,1,1,kb) ) .* art_3d ) ;
qf(1,:,:) =0;  qf(im,:,:)=0;
qf(:,1,:) =0;  qf(:,jm,:)=0;   
qf(:,:,1) =0;  qf(:,:,kb)=0;
return
