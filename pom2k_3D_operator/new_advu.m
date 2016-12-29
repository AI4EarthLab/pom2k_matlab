function [uf] = new_advu(advx,cor,dt,e_atmos,drhox,h,ub,u,v,w,egf,egb,etf,etb)
%load('grid.mat');load('para.mat');load('operator.mat');

global im jm kb aru_3d dti2 grav dy_3d dz_3d;

h_3d = repmat(h,1,1,kb); 
dt_3d = repmat(dt,1,1,kb); 
etb_3d = repmat(etb,1,1,kb);
cor_3d = repmat(cor,1,1,kb);
egf_3d = repmat(egf,1,1,kb);
egb_3d = repmat(egb,1,1,kb); 
e_atmos_3d = repmat(e_atmos,1,1,kb);
etf_3d = repmat(etf,1,1,kb);

uf=zeros(im,jm,kb);
uf = AXB(w) .* AZB(u);
bond = uf(im,:,:);
uf=DIVISION( AXB(etb_3d+h_3d).*aru_3d.*ub...
                -dti2*( advx+drhox-aru_3d.*AXB( cor_3d.*dt_3d.*AYF(v) )...
                +grav*AXB(dt_3d).*( DXB(egf_3d+egb_3d)+DXB(e_atmos_3d)*2.0 ).*AXB(dy_3d)/2.e0...
                -DZF(uf).*aru_3d./dz_3d )...
                ,AXB(etf_3d+h_3d).*aru_3d );
uf(im,:,:) = bond;

