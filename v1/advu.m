function [uf] = advu(advx,cor,dt_3d,e_atmos,drhox,ub,u,v,w,egf,egb,etf,etb)
global im kb dti2 grav h_3d;

uf=DIVISION( (AXB(repmat(etb,1,1,kb)+h_3d).*ub -dti2*( advx + drhox - AXB( repmat(cor,1,1,kb) .*dt_3d.*AYF(v) )...
   +grav*AXB(dt_3d).*( DXB( repmat(egf+egb,1,1,kb) )+DXB( repmat(e_atmos,1,1,kb) )*2.0 )/2.e0-DZF(AXB(w) .* AZB(u)))), ...
   AXB( repmat(etf,1,1,kb)+h_3d ) );
uf(:,:,kb)=0.e0;  bond=AXB(w) .* AZB(u); uf(im,:,:) = bond(im,:,:) ;   %add by hx

