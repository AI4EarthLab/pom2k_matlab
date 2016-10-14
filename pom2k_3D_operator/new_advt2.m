
function [fb,f,fclim,ff,sw] = advt2(fb,f,fclim,ff,xflux,yflux,etb,etf,w,dt,aam,h,u,v)
load('grid.mat');load('para.mat');load('operator.mat');
xmassflux=zeros(im,jm,kb);ymassflux=zeros(im,jm,kb);
xflux=zeros(im,jm,kb);yflux=zeros(im,jm,kb);zflux=zeros(im,jm,kb);
%
for k=1:kbm1
    xmassflux(:,:,k)=AXB1_XY(dy).*AXB1_XY(dt).*u(:,:,k);
    ymassflux(:,:,k)=AYB1_XY(dx).*AYB1_XY(dt).*v(:,:,k);
end
fb(:,:,kb)=fb(:,:,kbm1);
eta=etb;zwflux=w;fbmem=fb;
%
%     Start Smolarkiewicz scheme:
%
for itera=1:nitera
    %     Upwind advection scheme:
    tempx1=xmassflux+abs(xmassflux);tempx2=xmassflux-abs(xmassflux);
    tempy1=ymassflux+abs(ymassflux);tempy2=ymassflux-abs(ymassflux);
    tempz1=zwflux+abs(zwflux);tempz2=zwflux-abs(zwflux);
    for i=2:im
        xflux(i,:,:)=0.5e0*( tempx1(i,:,:).*fbmem(i-1,:,:)+tempx2(i,:,:).*fbmem(i,:,:) );
    end
    for j=2:jm
        yflux(:,j,:)=0.5e0*( tempy1(:,j,:).*fbmem(:,j-1,:)+tempy2(:,j,:).*fbmem(:,j,:) );
    end
%
    zflux(:,:,1)=0.e0;
    if(itera==1)
       zflux(:,:,1)=w(:,:,1).*f(:,:,1).*art;
    end
    zflux(:,:,kb)=0.e0;

    %
    for k=2:kbm1
        zflux(:,:,k)=0.5e0*( tempz1(:,:,k).*fbmem(:,:,k)+tempz2(:,:,k).*fbmem(:,:,k-1) ).*art;
    end
    %
    %     Add net advective fluxes and step forward in time:
    %
    temp=zeros(im,jm,kb);
    for j=2:jmm1 
    temp(:,j,:)=-DZF1_XZ( permute(zflux(:,j,:),[1,3,2]) );
    end
    bond=ff(im,:,:);
    for k=1:kbm1
    ff(:,:,k)=(fbmem(:,:,k).*(OP_L_XY*(h+eta)*OP_R_XY).*art...
              -dti2*( DXF2_XY(xflux(:,:,k))+DYF2_XY(yflux(:,:,k))+temp(:,:,k)/dz(k) ) )...
              ./((h+etf).*art);
    end
    ff(im,:,:)=bond;
    %
    %     %alculate antidiffusion velocity:
    %
   [xmassflux,ymassflux,zwflux,ff] = smol_adif(xmassflux,ymassflux,zwflux,ff,dt);
    eta=etf;
    %     End of Smolarkiewicz scheme
end

fb=fb-fclim;
%     Add net horizontal fluxes and step forward in time: 
%
for k=1:kbm1
    ff(:,:,k)=ff(:,:,k)...
              -dti2*( DXF2_XY( (-AXB1_XY(aam(:,:,k)).*AXB1_XY(h)*tprni.*DXB1_XY(fb(:,:,k)).*dum.*(AXB1_XY(dy)).*DIVISION(1.e0,AXB1_XY(dx))) )     ...
              +DYF2_XY( (-AYB1_XY(aam(:,:,k)).*AYB1_XY(h)*tprni.*DYB1_XY(fb(:,:,k)).*dvm.*(AYB1_XY(dx)).*DIVISION(1.e0,AYB1_XY(dy))) ) )...
              ./((h+etf).*art);
end
fb=fb+fclim;
return;
