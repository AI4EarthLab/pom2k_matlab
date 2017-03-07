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
global dx_3d dy_3d

curv= ( AYF(v) .* DXB(AXF(dy_3d)) - AXF(u) .* DYB(AYF(dx_3d)) ) ./ (dx_3d .* dy_3d);

advx  = DXB( AXF( AXB(dt_3d) .* u ) .* AXF(u) - dt_3d.*aam*2.0.*DXF(ub) ) + ...
        DYF( AXB( AYB(dt_3d) .* v ) .* AYB(u)  - AYB(AXB(dt_3d)) .* AYB(AXB(aam)).*( DYB(ub) + DXB(vb) ) ) - ...
        AXB(  curv .* dt_3d .* AYF(v) ); 
                 
advy  = DXF( AYB( AXB(dt_3d) .* u ) .* AXB(v) -AYB(AXB(dt_3d)) .* AYB(AXB(aam)).* (DYB(ub) + DXB(vb) ) ) + ...
        DYB( AYF( AYB(dt_3d) .* v ) .* AYF(v) - dt_3d.*aam*2.0.*DYF(vb) ) + ... 
        AYB(  curv .* dt_3d .* AXF(u) ); 
end
