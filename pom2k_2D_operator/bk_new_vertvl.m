function [xflux,yflux,w]=new_vertvl (w,dt,u,v,vfluxb,vfluxf,etf,etb,dti2)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  calculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('xyz.mat'); load('operator.mat');
xflux = zeros(im,jm,kb);
yflux = zeros(im,jm,kb);
temp=zeros(im,jm,kb);
for k=1:kbm1
    xflux(:,:,k) = AXB1_XY(dy) .* AXB1_XY(dt) .* u(:,:,k);
    yflux(:,:,k) = AYB1_XY(dx) .* AYB1_XY(dt) .* v(:,:,k); 
    %etf-etb在边界处不能有值
    temp(:,:,k)=dz(k).* ( ( DXF2_XY(xflux(:,:,k))+DYF2_XY(yflux(:,:,k)) )...
                ./(dx.*dy)+OP_L_XY*(etf-etb)/dti2 );
end
%
%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
%
w = zeros(im,jm,kb);
w(:,:,1)=0.5*(vfluxb+vfluxf);
for j=2:jmm1            
    w(:,j,:)=OP_L_XZ*repmat(w(:,j,1),1,kb)+SUMZ2( permute(temp(:,j,:),[1,3,2]) );
end
return



