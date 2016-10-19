function [vf] = advv(advy,cor,dt,e_atmos,drhoy,h,vb,vf,u,v,w,egf,egb,etf,etb)
load('grid.mat');load('para.mat');load('operator.mat');

h_3d = repmat(h,1,1,kb); 
dt_3d = repmat(dt,1,1,kb); 
dx_3d = repmat(dx,1,1,kb);
etb_3d = repmat(etb,1,1,kb);
arv_3d = repmat(arv,1,1,kb); 
cor_3d = repmat(cor,1,1,kb);
egf_3d = repmat(egf,1,1,kb);
egb_3d = repmat(egb,1,1,kb); 
etf_3d = repmat(etf,1,1,kb);
e_atmos_3d = repmat(e_atmos,1,1,kb);


vf=zeros(im,jm,kb);
vf = AYB1(w) .* AZB2(v);
vf = DIVISION( AYB2(etb_3d+h_3d).*arv_3d.*vb...
              -dti2*( advy+drhoy+arv_3d.*AYB2( cor_3d.*dt_3d.*AXF2(u) )...
              +grav*AYB2(dt_3d).*( DYB2(egf_3d+egb_3d)+DYB2(e_atmos_3d)*2.0 ).*AYB2(dx_3d)/2.e0...
              -DZF1(vf).*arv_3d./dz_3d )...
              ,AYB2(etf_3d+h_3d).*arv_3d );
