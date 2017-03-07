function [vf] = advv(advy,cor,dt_3d,e_atmos,drhoy,vb,u,v,w,egf,egb,etf,etb)
global kb dti2 grav h_3d

vf = ( AYB( repmat(etb,1,1,kb) +h_3d).*vb-dti2*( advy+drhoy+AYB( repmat(cor,1,1,kb) .*dt_3d.*AXF(u) )...
     +grav*AYB(dt_3d).*( DYB( repmat(egf+egb,1,1,kb) )+DYB( repmat(e_atmos,1,1,kb) )*2.0 )/2.e0-DZF(AYB(w) .* AZB(v)))) ...
     ./AYB( repmat(etf,1,1,kb) +h_3d);
vf(:,:,kb)=0.e0;    %add by hx 
return
