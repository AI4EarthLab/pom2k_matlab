function [qb,q,qf,xflux,yflux]=advq(qb,q,qf,xflux,yflux,...
    dt,dx,dy,dz,u,v,w,aam,h,dum,dvm,art,etb,etf,im,jm,imm1,jmm1,kbm1,dti2)

% **********************************************************************
% *                                                                    *
% * FUN%TION    :  %alculates horizontal advection and diffusion, and  *
% *                vertical advection for turbulent quantities.        *
% *                                                                    *
% **********************************************************************
%
%     Do horizontal advection:
%

%%%%%----------------------------------------------------------------------
xflux = zeros(im,jm,kb);
yflux = zeros(im,jm,kb);

%%%%%%% which method is better
for j=2:jm
    xflux(:,j,:)=AXB1_XZ(permute(q(:,j,:),[1,3,2])) .* AXB1_XZ(repmat(dt(:,j),1,kb)) .* AZB2_XZ(permute(u(:,j,:),[1,3,2]));
end

for i=2:im
    yflux(i,:,:)=AYB1_YZ(permute(q(i,:,:),[2,3,1])) .* AYB1_YZ(repmat(dt(i,:)',1,kb)) .* AZB2_YZ(permute(v(i,:,:),[2,3,1]));
end

%
%     Do horizontal diffusion:
%
for k=2:kbm1
    xflux(:,:,k) = xflux(:,:,k)...
        -0.5*( AXB1_XY(aam(:,:,k))+AXB1_XY(aam(:,:,k-1)) )...
        .* AXB1_XY(h) .* DXB1_XY(qb(:,:,k)).*dum ./ AXB1_XY(dx);
    yflux(:,:,k) = yflux(:,:,k)...
        -0.5*( AYB1_XY(aam(:,:,k))+AYB1_XY(aam(:,:,k-1)) )...
        .* AYB1_XY(h) .* DYB1_XY(qb(:,:,k)).*dvm ./ AYB1_XY(dy);
    xflux(:,:,k)=AXB1_XY(dy) .* xflux(:,:,k);
    yflux(:,:,k)=AYB1_XY(dx) .* yflux(:,:,k);
end
xflux(find(isnan(xflux)==1)) = 0;
yflux(find(isnan(yflux)==1)) = 0;

%
%     do vertical advection: add flux terms, then step forward in time:
%
qf=zeros(im,jm,kb);

% !!!pay attention to the boundary
for k=2:kbm1
    qf(:,:,k)=OP_L_XY*(w(:,:,k-1).*q(:,:,k-1)-w(:,:,k+1).*q(:,:,k+1))*OP_R_XY  ...
                .*art/(dz(k)+dz(k-1))     ...
                +DXF2_XY(xflux(:,:,k))+DYF2_XY(yflux(:,:,k));   
     qf(:,:,k)=OP_L_XY*( (h+etb) .* art .* qb(:,:,k)-dti2.*qf(:,:,k) )*OP_R_XY ...
         ./((h+etf).*art);  
end

return
end