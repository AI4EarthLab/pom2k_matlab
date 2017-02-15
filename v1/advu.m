function [uf] = advu(advx,cor,dt_3d,e_atmos,drhox,h_3d,ub,u,v,w,egf,egb,etf,etb)
global im jm kb dti2 grav gs;

etb_3d = repmat(etb,1,1,kb);
cor_3d = repmat(cor,1,1,kb);
egf_3d = repmat(egf,1,1,kb);
egb_3d = repmat(egb,1,1,kb);
etf_3d = repmat(etf,1,1,kb);
e_atmos_3d = repmat(e_atmos,1,1,kb);

uf=create_field(zeros(im,jm,kb),gs,2);
uf = AXB(w) .* AZB(u);
bond = uf(im,:,:);
uf=(AXB(etb_3d+h_3d).*ub-dti2*( advx + drhox - AXB( cor_3d.*dt_3d.*AYF(v) )...
   +grav*AXB(dt_3d).*( DXB(egf_3d+egb_3d)+DXB(e_atmos_3d)*2.0 )/2.e0-DZF(uf)))./AXB(etf_3d+h_3d);
uf(im,:,:) = bond;

