function [qf]=advq(qb,q,dt_3d,u,v,w,aam,etb,etf)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates horizontal advection and diffusion, and  *
% *                vertical advection for turbulent quantities.        *
% *                                                                    *
% **********************************************************************
global kb dum_3d dvm_3d dti2 h_3d

qf= ( (h_3d+repmat(etb,1,1,kb)) .* qb   ...
      -dti2*( -DZB(AZF(w.*q)) + DXF(AXB(q) .* AXB(dt_3d) .* AZB(u) -AZB( AXB( aam )).*AXB(h_3d) .*DXB( qb ).* dum_3d)  ...
      +DYF(AYB(q) .* AYB(dt_3d) .* AZB(v) -AZB( AYB( aam )).*AYB(h_3d) .*DYB( qb ).* dvm_3d)) )./  ( h_3d+repmat(etf,1,1,kb) ) ;
return
