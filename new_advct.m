function [advx,advy]= advct(u,v,dx,dy,dt,aam,ub,vb,aru,arv,im,jm,kb,imm1,jmm1,kbm1)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates the horizontal portions of momentum      *
% *                advection well in advance of their use in advu and  *
% *                advv so that their vertical integrals (created in   *
% *                the main program) may be used in the external (2-D) *
% *                mode calculation.                                   *
% *                                                                    *
% **********************************************************************
%
%
curv=zeros(im,jm,kb);
advx=zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);

[AXF, AXB, AYF, AYB, DXF, DXB, DYF, DYB]=operator(im,jm);
DXBAXF= DXB * AXF;
DYBAYF= AYF * DYB;   %\delta_y avg(F)^x

for k=1:kbm1
    curv(:,:,k)=(v(:,:,k)*AYF) .* (DXBAXF *dy)  - (AXF*u(:,:,k)) .* (dx* DYBAYF) ./ (dx .* dy);
end

%     Calculate x-component of velocity advection:
%
%     Calculate horizontal advective fluxes:
for k=1:kbm1
    xflux(:,:,k) = ( AXF*( (AXB*dt) .* u(:,:,k) ) ) .* ( AXF*u(:,:,k) );
end
xflux(1, :, :)=0;
xflux(im,:, :)=0;

for k=1:kbm1
    yflux(:,:,k)=( AXB*(dt * AYB) .* v(:,:,k) ) .* ( u(:,:,k) * AYB ) ;
end
yflux(1,:,k)=0;
yflux(:,1,k)=0;


%    Add horizontal diffusive fluxes:
dtaam=zeros(im,jm);
for k=1:kbm1
    xflux(:,:,k)=xflux(:,:,k) - dt.*aam(:,:,k)*2.0.*(DXF*ub(:,:,k))./dx;
    dtaam(:,:)=(AXB*dt*AYB).*(AXB*aam(:,:,k)*AYB);
    yflux(:,:,k)=yflux(:,:,k)-dtaam(:,:).*((ub(:,:,k)*DYB)./(AXB*dy*AYB)+(DXB*vb(:,:,k))./(AXB*dx*AYB));
    xflux(:,:,k)=dy.*xflux(:,:,k);
    yflux(:,:,k)=(AXB*dx*AYB).*yflux(:,:,k);
end
xflux(1, :, :)=0;
xflux(im,:, :)=0;
xflux(:, 1, :)=0;

yflux(1, :, :)=0;
yflux(im,:, :)=0;
yflux(:, 1, :)=0;

%     for horizontal advection:
for k=1:kbm1
    advx(:,:,k)=DXB*xflux(:,:,k)+ yflux(:,:,k)*DYF ;
end
advx(1,:,:)=0;
advx(im,:,:)=0;
advx(:,1,:)=0;
advx(:,jm,:)=0;

Work=zeros(im,jm);
for k=1:kbm1
   Work(:,:) = aru.*( AXB * (curv(:,:,k) .* dt .* (v(:,:,k)*AYF)) );
   advx(3:imm1, 2:jmm1, k) = advx(3:imm1, 2:jmm1, k)-Work(3:imm1, 2:jmm1);
end

%
%-----------------------------------------------------------------------
%
advy=zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);
%
%     Calculate y-component of velocity advection:
%
%     Calculate horizontal advective fluxes:
for k=1:kbm1
    xflux(:,:,k) = ( ( (AXB*dt) .* u(:,:,k) ) * AYB ) .* ( AXB*v(:,:,k) );
end
xflux(1, :, :)=0;
xflux(:, 1, :)=0;

for k=1:kbm1
    yflux(:,:,k)=( (dt * AYB) .* v(:,:,k) * AYF ) .* ( v(:,:,k) * AYF ) ;
end
yflux(:,1,:)=0;
 
%    Add horizontal diffusive fluxes:
dtaam=zeros(im,jm);
for k=1:kbm1
    dtaam(:,:)=(AXB*dt*AYB).*(AXB*aam(:,:,k)*AYB);
    xflux(:,:,k)=xflux(:,:,k)-dtaam(:,:).*((ub(:,:,k)*DYB)./(AXB*dy*AYB)+(DXB*vb(:,:,k))./(AXB*dx*AYB));
    yflux(:,:,k)=yflux(:,:,k) - dt.*aam(:,:,k)*2.0.*(DXF*ub(:,:,k))./dy;
    xflux(:,:,k)=(AXB*dy*AYB).*xflux(:,:,k);
    yflux(:,:,k)=dx.*yflux(:,:,k);
end
xflux(1, :, :)=0;
xflux(:, 1, :)=0;
xflux(:,jm, :)=0;
yflux(1, :, :)=0;
yflux(:, 1, :)=0;
yflux(:,jm, :)=0;

%     for horizontal advection:
for k=1:kbm1
    advy(:,:,k)=DXF*xflux(:,:,k)+ yflux(:,:,k)*DYB ;
end
advy(1,:,:)=0;
advy(im,:,:)=0;
advy(:,1,:)=0;
advy(:,jm,:)=0;

Work=zeros(im,jm);
for k=1:kbm1
   Work(:,:) = arv.*( ( curv(:,:,k) .* dt .* (AXF*u(:,:,k)) ) * AYB );
   advy(2:imm1, 3:jmm1, k) = advy(2:imm1, 3:jmm1, k)-Work(2:imm1, 3:jmm1);
end