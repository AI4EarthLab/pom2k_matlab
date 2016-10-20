function [uf] = advu(advx,cor,dt,e_atmos,drhox,h,ub,uf,u,v,w,egf,egb,etf,etb)
load('grid.mat');load('para.mat');load('operator.mat');

h_3d = repmat(h,1,1,kb); 
dt_3d = repmat(dt,1,1,kb); 
dy_3d = repmat(dy,1,1,kb);
etb_3d = repmat(etb,1,1,kb);
aru_3d = repmat(aru,1,1,kb); 
cor_3d = repmat(cor,1,1,kb);
egf_3d = repmat(egf,1,1,kb);
egb_3d = repmat(egb,1,1,kb); 
e_atmos_3d = repmat(e_atmos,1,1,kb);
etf_3d = repmat(etf,1,1,kb);

uf=zeros(im,jm,kb);
uf = AXB1(w) .* AZB2(u);
bond = uf(im,:,:);
uf=DIVISION( AXB2(etb_3d+h_3d).*aru_3d.*ub...
                -dti2*( advx+drhox-aru_3d.*AXB2( cor_3d.*dt_3d.*AYF2(v) )...
                +grav*AXB2(dt_3d).*( DXB2(egf_3d+egb_3d)+DXB2(e_atmos_3d)*2.0 ).*AXB2(dy_3d)/2.e0...
                -DZF1(uf).*aru_3d./dz_3d )...
                ,AXB2(etf_3d+h_3d).*aru_3d );
uf(im,:,:) = bond;

