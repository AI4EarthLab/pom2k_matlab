function [advx,advy]=advct(u,v,dt_3d,aam,ub,vb)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates the horizontal portions of momentum      *
% *                advection well in advance of their use in advu and  *
% *                advv so that their vertical integrals (created in   *
% *                the main program) may be used in the external (2-D) *
% *                mode calculation.                                   *
% *                                                                    *
% **********************************************************************
global im jm kb dx_3d dy_3d gs;

curv    = create_field(zeros(im,jm,kb),    gs, 3);
xflux   = create_field(zeros(im,jm,kb),    gs, 3);
yflux   = create_field(zeros(im,jm,kb),    gs, 3);

%curv= (AYF(v) .* DXB(AXF(dy_3d)) - AXF(u) .* DYB(AYF(dx_3d)) ) ./ (dx_3d .* dy_3d);
curv= AYF(v) .* DXB(AXF(dy_3d))./dy_3d -  AXF(u) .* DYB(AYF(dx_3d))./dx_3d;

% Calculate x-component of velocity advection:
%xflux = dy_3d .* (AXF( AXB(dt_3d) .* u ) .* AXF(u) - dt_3d.*aam*2.0.*DXF(ub)./dx_3d );                   
%yflux = AYB(AXB(dx_3d)) .* ( ( AXB( AYB(dt_3d) .* v ) .* AYB(u) ) - AYB(AXB(dt_3d)) .* AYB(AXB(aam)) ...
%                 .*( DIVISION( DYB(ub), AYB(AXB(dy_3d)) ) + DIVISION( DXB(vb), AYB(AXB(dx_3d)) ) ) ) ;
%advx = DXB( xflux ) + DYF( yflux ) - aru_3d .* AXB(  curv .* dt_3d .* AYF(v) ); 
xflux = AXF( AXB(dt_3d) .* u ) .* AXF(u) - dt_3d.*aam*2.0.*DXF(ub) ;
yflux = AXB( AYB(dt_3d) .* v ) .* AYB(u)  - AYB(AXB(dt_3d)) .* AYB(AXB(aam)).*( DYB(ub) + DXB(vb) ) ;
advx  = DXB( xflux ) + DYF( yflux ) - AXB(  curv .* dt_3d .* AYF(v) ); 

% Calculate y-component of velocity advection:
% % % xflux = AYB(AXB(dy_3d)) .* ( ( AYB( AXB(dt_3d) .* u ) .* AXB(v) ) -AYB(AXB(dt_3d)) .* AYB(AXB(aam)) ...
% % %                  .* ( DIVISION( DYB(ub), AYB(AXB(dy_3d)) ) + DIVISION( DXB(vb), AYB(AXB(dx_3d)) ) ) ) ;
% % % yflux = dx_3d .* (AYF( AYB(dt_3d) .* v ) .* AYF(v) - dt_3d.*aam*2.0.*DYF(vb)./dy_3d );                   
% % % advy=DXF( xflux ) + DYB( yflux ) + arv_3d .* AYB(  curv .* dt_3d .* AXF(u) ); 

xflux = AYB( AXB(dt_3d) .* u ) .* AXB(v) -AYB(AXB(dt_3d)) .* AYB(AXB(aam)).* (DYB(ub) + DXB(vb) )  ;
yflux = AYF( AYB(dt_3d) .* v ) .* AYF(v) - dt_3d.*aam*2.0.*DYF(vb);                   
advy  = DXF( xflux ) + DYB( yflux ) + AYB(  curv .* dt_3d .* AXF(u) ); 

end
