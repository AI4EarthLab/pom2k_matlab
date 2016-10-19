function [fb,f,fclim,ff,xflux,yflux,zflux]=advt1(fb,f,fclim,ff,xflux,yflux,...
                                                zflux,dt,u,v,aam,tprni,dum,dvm,w,art,etb,etf,dti2,dx,dy,dz,h,im,jm,kb,imm1,jmm1,kbm1)


% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is centred scheme, as originally provide in    *
% *                POM (previously called advt).                       *
% *                                                                    *
% **********************************************************************
%
%
%

for j=1:jm
    for i=1:im
        f(i,j,kb)=f(i,j,kbm1);
        fb(i,j,kb)=fb(i,j,kbm1);
    end
end
%
%     for advective fluxes:
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=.25e0*((dt(i,j)+dt(i-1,j)) ...
                *(f(i,j,k)+f(i-1,j,k))*u(i,j,k));
            yflux(i,j,k)=.25e0*((dt(i,j)+dt(i,j-1))     ...
                *(f(i,j,k)+f(i,j-1,k))*v(i,j,k));
        end
    end
end
%
%     Add diffusive fluxes:
%

fb = fb-fclim;
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=xflux(i,j,k)     ...
                -.5e0*(aam(i,j,k)+aam(i-1,j,k))     ...
                *(h(i,j)+h(i-1,j))*tprni     ...
                *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)     ...
                /(dx(i,j)+dx(i-1,j));
            yflux(i,j,k)=yflux(i,j,k)     ...
                -.5e0*(aam(i,j,k)+aam(i,j-1,k))     ...
                *(h(i,j)+h(i,j-1))*tprni     ...
                *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)     ...
                /(dy(i,j)+dy(i,j-1));
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k);
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k);
        end
    end
end
%
fb = fb+fclim;
%
%     for vertical advection:
%
for j=2:jmm1
    for i=2:imm1
        zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j);
        zflux(i,j,kb)=0.e0;
    end
end
%
for k=2:kbm1
    for j=2:jmm1
        for i=2:imm1
            zflux(i,j,k)=.5e0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j);
        end
    end
end
%
%     Add net horizontal fluxes and then step forward in time:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)     ...
                +yflux(i,j+1,k)-yflux(i,j,k)     ...
                +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k);
            %
            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)     ...
                -dti2*ff(i,j,k))     ...
                /((h(i,j)+etf(i,j))*art(i,j));
        end
    end
end
return 
end

