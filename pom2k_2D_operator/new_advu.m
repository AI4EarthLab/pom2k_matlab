function [uf] = advu(advx,cor,dt,e_atmos,drhox,h,ub,uf,u,v,w,egf,egb,etf,etb)
load('grid.mat');load('depth.mat');load('para.mat');load('operator.mat');load('xyz.mat');load('masks.mat');
uf=zeros(im,jm,kb);
%
for j=1:jm
    uf(:,j,:)=AXB1_XZ( permute(w(:,j,:),[1,3,2]) ) .* AZB2_XZ( permute(u(:,j,:),[1,3,2]) );
end

for k=1:kbm1
    uf(:,:,k)=( AXB2_XY(etb+h).*aru.*ub(:,:,k)...
                -dti2*( advx(:,:,k)+drhox(:,:,k)-aru.*AXB2_XY( cor.*dt.*AYF2_XY( v(:,:,k) ) )...
                +grav*AXB2_XY(dt).*( DXB2_XY(egf+egb)+DXB2_XY(e_atmos)*2.0 ).*AXB2_XY(dy)/2.e0...
                +(uf(:,:,k)-uf(:,:,k+1)).*aru/dz(k) ) )...
                .*DIVISION( 1.e0,AXB2_XY(etf+h).*aru );
end
return