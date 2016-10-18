function [advx,advy]=new_advct(u,v,dt_3d,aam,ub,vb)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates the horizontal portions of momentum      *
% *                advection well in advance of their use in advu and  *
% *                advv so that their vertical integrals (created in   *
% *                the main program) may be used in the external (2-D) *
% *                mode calculation.                                   *
% *                                                                    *
% **********************************************************************
load('grid.mat');
curv=zeros(im,jm,kb);  advx=zeros(im,jm,kb); advy=zeros(im,jm,kb);
xflux=zeros(im,jm,kb); yflux=zeros(im,jm,kb);

curv= (AYF1(v) .* DXB2(AXF1(dy_3d)) - AXF1(u) .* DYB2(AYF1(dx_3d)) ) ./ (dx_3d .* dy_3d);
% Calculate x-component of velocity advection:
xflux = dy_3d .* (AXF1( AXB1(dt_3d) .* u ) .* AXF2(u) - dt_3d.*aam*2.0.*DXF2(ub)./dx_3d );                   
yflux = AYB1(AXB1(dx_3d)) .* ( ( AXB1( AYB1(dt_3d) .* v ) .* AYB1(u) ) - AYB1(AXB1(dt_3d)) .* AYB1(AXB1(aam)) ...
                 .*( DIVISION( DYB1(ub), AYB1(AXB1(dy_3d)) ) + DIVISION( DXB1(vb), AYB1(AXB1(dx_3d)) ) ) ) ;
advx=DXB2( xflux ) + DYF2( yflux ) - aru_3d .* AXB2(  curv .* dt_3d .* AYF2(v) ); 
% Calculate y-component of velocity advection:
xflux = AYB1(AXB1(dy_3d)) .* ( ( AYB1( AXB1(dt_3d) .* u ) .* AXB1(v) ) -AYB1(AXB1(dt_3d)) .* AYB1(AXB1(aam)) ...
                 .* ( DIVISION( DYB1(ub), AYB1(AXB1(dy_3d)) ) + DIVISION( DXB1(vb), AYB1(AXB1(dx_3d)) ) ) ) ;
yflux = dx_3d .* (AYF1( AYB1(dt_3d) .* v ) .* AYF2(v) - dt_3d.*aam*2.0.*DYF2(vb)./dy_3d );                   
advy=DXF2( xflux ) + DYB2( yflux ) - aru_3d .* AYB2(  curv .* dt_3d .* AXF2(u) ); 

end