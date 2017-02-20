function [vf] = advv(advy,cor,dt_3d,e_atmos,drhoy,h_3d,vb,u,v,w,egf,egb,etf,etb)
global im jm kb dti2 grav gs

etb_3d = repmat(etb,1,1,kb);
cor_3d = repmat(cor,1,1,kb);
egf_3d = repmat(egf,1,1,kb);
egb_3d = repmat(egb,1,1,kb); 
etf_3d = repmat(etf,1,1,kb);
e_atmos_3d = repmat(e_atmos,1,1,kb);

vf=create_field(zeros(im,jm,kb),gs,1);
vf = AYB(w) .* AZB(v);
% vf = DIVISION( AYB(etb_3d+h_3d).*arv_3d.*vb...
%               -dti2*( advy+drhoy+arv_3d.*AYB( cor_3d.*dt_3d.*AXF(u) )...
%               +grav*AYB(dt_3d).*( DYB(egf_3d+egb_3d)+DYB(e_atmos_3d)*2.0 ).*AYB(dx_3d)/2.e0...
%               -DZF(vf).*arv_3d./dz_3d )...
%               ,AYB(etf_3d+h_3d).*arv_3d );
vf = ( AYB(etb_3d+h_3d).*vb-dti2*( advy+drhoy+AYB( cor_3d.*dt_3d.*AXF(u) )...
     +grav*AYB(dt_3d).*( DYB(egf_3d+egb_3d)+DYB(e_atmos_3d)*2.0 )/2.e0-DZF(vf)))./AYB(etf_3d+h_3d);
vf(:,:,kb)=0.e0;    %add by hx 
return
