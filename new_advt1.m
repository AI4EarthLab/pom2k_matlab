function [fb,f,fclim,ff]=advt1(fb,f,fclim,ff,u,v,aam,w,h)
load('grid.mat');load('depth.mat');load('para.mat');load('operator.mat');load('xyz.mat');load('masks.mat');
xflux=zeros(im,jm,kb);yflux=zeros(im,jm,kb);
f(:,:,kb)=f(:,:,kbm1);fb(:,:,kb)=fb(:,:,kbm1);
zflux(:,:,1)=f(:,:,1).*w(:,:,1).*art;     
zflux(:,:,kb)=0.e0;
% % %如果没有上述边界条件，以下两个循环可以合并成一个
for j=2:jmm1
    zflux(:,j,:)=AZB2_XZ( permute(f(:,j,:),[1,3,2]) ).*permute(w(:,j,:),[1,3,2]).*repmat(art(:,j),1,kb);
end
temp=zeros(im,jm,kb);
for j=2:jmm1 
    temp(:,j,:)=-DZF1_XZ( permute(zflux(:,j,:),[1,3,2]) );
end
%%%%%%
bond=ff(im,:,:);
for k=1:kbm1
    ff(:,:,k)=( fb(:,:,k).*( OP_L_XY*(h+etb)*OP_R_XY ).*art...
              -dti2* ( DXF2_XY( AXB1_XY(dy).*( AXB1_XY(dt).*AXB1_XY(f(:,:,k)).*u(:,:,k)-AXB1_XY(aam(:,:,k)).*AXB1_XY(h)*tprni.*DXB1_XY(fb(:,:,k)-fclim(:,:,k)).*dum.*DIVISION( 1.0,AXB1_XY(dx) ) ) )...
              +DYF2_XY( AYB1_XY(dx).*( AYB1_XY(dt).*AYB1_XY(f(:,:,k)).*v(:,:,k)-AYB1_XY(aam(:,:,k)).*AYB1_XY(h)*tprni.*DYB1_XY(fb(:,:,k)-fclim(:,:,k)).*dvm.*DIVISION( 1.0,AYB1_XY(dy) ) ) )...
              +temp(:,:,k)/dz(k) ) )...
              ./((h+etf).*art);
end
ff(im,:,:)=bond;
return