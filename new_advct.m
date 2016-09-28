function [advx,advy]=new_advct(u,v,dx,dy,dt,aam,ub,vb,aru,arv)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates the horizontal portions of momentum      *
% *                advection well in advance of their use in advu and  *
% *                advv so that their vertical integrals (created in   *
% *                the main program) may be used in the external (2-D) *
% *                mode calculation.                                   *
% *                                                                    *
% **********************************************************************
%
%
load('grid.mat');

curv=zeros(im,jm,kb);
advx=zeros(im,jm,kb);
advy=zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);

for k=1:kb
    curv(:,:,k)= (AYF1_XY(v(:,:,k)) .* DXB2_XY(AXF1_XY(dy)) - AXF1_XY(u(:,:,k)) .* DYB2_XY(AYF1_XY(dx)) ) ./ (dx .* dy);
%     Calculate x-component of velocity advection:
    xflux(:,:,k) = dy .* (AXF1_XY( AXB1_XY(dt) .* u(:,:,k) ) .* AXF2_XY(u(:,:,k)) ...
                    - dt.*aam(:,:,k)*2.0.*DXF2_XY(ub(:,:,k))./dx );                   
    yflux(:,:,k) = AYB1_XY(AXB1_XY(dx)) .* ( ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k)) ) ...
                     -AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
                     .*( DIVISION( DYB1_XY(ub(:,:,k)), AYB1_XY(AXB1_XY(dy)) ) ...
                       + DIVISION( DXB1_XY(vb(:,:,k)), AYB1_XY(AXB1_XY(dx)) ) ) ) ;
    advx(:,:,k)=DXB2_XY( xflux(:,:,k) ) + DYF2_XY( yflux(:,:,k) )...
                    - aru .* AXB2_XY(  curv(:,:,k) .* dt .* AYF2_XY(v(:,:,k)) ); 
%     Calculate y-component of velocity advection:
    xflux(:,:,k) = AYB1_XY(AXB1_XY(dy)) .* ( ( AYB1_XY( AXB1_XY(dt) .* u(:,:,k) ) .* AXB1_XY(v(:,:,k)) ) ...
                     -AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
                     .*( DIVISION( DYB1_XY(ub(:,:,k)), AYB1_XY(AXB1_XY(dy)) ) ...
                       + DIVISION( DXB1_XY(vb(:,:,k)), AYB1_XY(AXB1_XY(dx)) ) ) ) ;
    yflux(:,:,k) = dx .* (AYF1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYF2_XY(v(:,:,k)) ...
                    - dt.*aam(:,:,k)*2.0.*DYF2_XY(vb(:,:,k))./dy );                   
    advy(:,:,k)=DXF2_XY( xflux(:,:,k) ) + DYB2_XY( yflux(:,:,k) )...
                    - aru .* AYB2_XY(  curv(:,:,k) .* dt .* AXF2_XY(u(:,:,k)) ); 
end