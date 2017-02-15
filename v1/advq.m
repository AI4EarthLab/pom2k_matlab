function [qf]=advq(qb,q,dt_3d,u,v,w,aam,h,etb,etf,dti2)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates horizontal advection and diffusion, and  *
% *                vertical advection for turbulent quantities.        *
% *                                                                    *
% **********************************************************************
%load('grid.mat');load('operator.mat');
global im jm kb dum_3d dvm_3d gs

xflux = Field(zeros(im,jm,kb),gs,6);
yflux = Field(zeros(im,jm,kb),gs,5);
qf=Field(zeros(im,jm,kb),gs,7);

h_3d =repmat(h,1,1,kb);

%     Do horizontal advection and diffusion:
%xflux=AXB(dy_3d) .* ( AXB(q) .* AXB(dt_3d) .* AZB(u) -AZB( AXB( aam ))...
%         .*AXB(h_3d) .*DXB( qb ) .* dum_3d .* DIVISION( 1.0,AXB( dx_3d ) ) );
xflux=( AXB(q) .* AXB(dt_3d) .* AZB(u) -AZB( AXB( aam )).*AXB(h_3d) .*DXB( qb ).* dum_3d );

%     do vertical advection: add flux terms, then step forward in time:
yflux=( AYB(q) .* AYB(dt_3d) .* AZB(v) -AZB( AYB( aam )).*AYB(h_3d) .*DYB( qb ).* dvm_3d);

qf= ( (h_3d+repmat(etb,1,1,kb)) .* qb   ...
      -dti2*( -DZB(AZF(w.*q)) + DXF(xflux)+DYF(yflux)) ...
    )./  ( h_3d+repmat(etf,1,1,kb) ) ;
qf(1,:,:) =0;  qf(im,:,:)=0;
qf(:,1,:) =0;  qf(:,jm,:)=0;   
qf(:,:,1) =0;  qf(:,:,kb)=0;
return
