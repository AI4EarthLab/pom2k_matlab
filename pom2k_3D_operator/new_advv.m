function [vf] = advv(advy,cor,dt,e_atmos,drhoy,h,vb,vf,u,v,w,egf,egb,etf,etb)
load('grid.mat');load('para.mat');load('operator.mat');

vf=zeros(im,jm,kb);
%
for i=1:im
    vf(i,:,:)=AYB1_YZ( permute(w(i,:,:),[2,3,1]) ) .* AZB2_YZ( permute(v(i,:,:),[2,3,1]) );
end
for k=1:kbm1
    vf(:,:,k)=( AYB2_XY(etb+h).*arv.*vb(:,:,k)...
              -dti2*( advy(:,:,k)+drhoy(:,:,k)+arv.*AYB2_XY( cor.*dt.*AXF2_XY( u(:,:,k) ) )...
              +grav*AYB2_XY(dt).*( DYB2_XY(egf+egb)+DYB2_XY(e_atmos)*2.0 ).*AYB2_XY(dx)/2.e0...
              +(vf(:,:,k)-vf(:,:,k+1)).*arv/dz(k) ) )...
              .*DIVISION( 1.e0,AYB2_XY(etf+h).*arv );
end
