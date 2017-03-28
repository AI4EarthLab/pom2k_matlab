function [fb,ff,xflux,yflux,zflux] = advt2(fb,f,fclim,ff,nitera,sw,etb,etf,w,art,dt,aam,tprni,h,dum,dvm,u,v,aru,arv)
global im jm kb imm1 jmm1 kbm1 dti2 dx dy dz;
%     calculate horizontal mass fluxes:
xmassflux=zeros(im,jm,kb);      ymassflux=zeros(im,jm,kb);      xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);          zflux=zeros(im,jm,kb);

for k=1:kbm1
    for j=2:jmm1
        for i=2:im
            xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j))*(dt(i-1,j)+dt(i,j))*u(i,j,k);
        end
    end
    %
    for j=2:jm
        for i=2:imm1
            ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j))*(dt(i,j-1)+dt(i,j))*v(i,j,k);
        end
    end
end
%
for j=1:jm
    for i=1:im
        fb(i,j,kb)=fb(i,j,kbm1);
    end
end
%
eta=etb;
zwflux=w;
fbmem=fb;
%
%     Start Smolarkiewicz scheme:
%
for itera=1:nitera
    %     Upwind advection scheme:

    for k=1:kbm1
        for j=2:jm
            for i=2:im
                xflux(i,j,k)=0.5e0 *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))     ...
                    *fbmem(i-1,j,k)+(xmassflux(i,j,k)-abs(xmassflux(i,j,k)))*fbmem(i,j,k));

                yflux(i,j,k)=0.5e0  *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))     ...
                    *fbmem(i,j-1,k)+ (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))*fbmem(i,j,k));
            end
        end
    end
    %
    for j=2:jmm1
        for i=2:imm1
            zflux(i,j,1)=0.e0;
            if(itera==1)
                zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j);
            end
            zflux(i,j,kb)=0.e0;
        end
    end
    %
    for k=2:kbm1
        for j=2:jmm1
            for i=2:imm1
                zflux(i,j,k)=0.5e0 *((zwflux(i,j,k)+abs(zwflux(i,j,k))) *fbmem(i,j,k)+     ...
                    (zwflux(i,j,k)-abs(zwflux(i,j,k)))*fbmem(i,j,k-1));
                zflux(i,j,k)=zflux(i,j,k)*art(i,j);
            end
        end
    end
    %
    %     Add net advective fluxes and step forward in time:
    %
    for k=1:kbm1
        for j=2:jmm1
            for i=2:imm1
                ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k)+(zflux(i,j,k)-zflux(i,j,k+1))/dz(k);
                ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)-dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j));
            end
        end
    end
    %
    %     %alculate antidiffusion velocity:
    %
  
   [xmassflux,ymassflux,zwflux,ff,sw] = smol_adif(xmassflux,ymassflux,zwflux,ff,sw,aru,arv,dt);
    %
    eta=etf;
    fbmem=ff;
    %     End of Smolarkiewicz scheme
end
%
%     Add horizontal diffusive fluxes:
%
fb=fb-fclim;
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k));
            ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k));
        end
    end
end
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni     ...
                *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)     ...
                *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j));
            yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni     ...
                *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)     ...
                *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1));
        end
    end
end

fb=fb+fclim;
%
%     Add net horizontal fluxes and step forward in time:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)     ...
                +yflux(i,j+1,k)-yflux(i,j,k)) /((h(i,j)+etf(i,j))*art(i,j));
        end
    end
end

end
