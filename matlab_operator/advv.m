function [vf] = advv(advy,cor,dt,e_atmos,drhoy,h,vb,u,v,w,egf,egb,etf,etb)
%load('grid.mat');load('para.mat');load('operator.mat');

global im jm kb arv_3d dti2 grav dx_3d dz_3d

h_3d = repmat(h,1,1,kb); 
dt_3d = repmat(dt,1,1,kb); 
etb_3d = repmat(etb,1,1,kb);
cor_3d = repmat(cor,1,1,kb);
egf_3d = repmat(egf,1,1,kb);
egb_3d = repmat(egb,1,1,kb); 
etf_3d = repmat(etf,1,1,kb);
e_atmos_3d = repmat(e_atmos,1,1,kb);


vf=zeros(im,jm,kb);
vf = AYB(w) .* AZB(v);
vf = DIVISION( AYB(etb_3d+h_3d).*arv_3d.*vb...
              -dti2*( advy+drhoy+arv_3d.*AYB( cor_3d.*dt_3d.*AXF(u) )...
              +grav*AYB(dt_3d).*( DYB(egf_3d+egb_3d)+DYB(e_atmos_3d)*2.0 ).*AYB(dx_3d)/2.e0...
              -DZF(vf).*arv_3d./dz_3d )...
              ,AYB(etf_3d+h_3d).*arv_3d );
              
return
