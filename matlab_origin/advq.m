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

for k=2:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=0.125*(q(i,j,k)+q(i-1,j,k))     ...
                *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1));
            yflux(i,j,k)=0.125*(q(i,j,k)+q(i,j-1,k))     ...
                *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1));
        end
    end
end
%
%     Do horizontal diffusion:
%
for k=2:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=xflux(i,j,k)     ...
                -0.25*(aam(i,j,k)+aam(i-1,j,k)     ...
                +aam(i,j,k-1)+aam(i-1,j,k-1))     ...
                *(h(i,j)+h(i-1,j))     ...
                *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)     ...
                /(dx(i,j)+dx(i-1,j));
            yflux(i,j,k)=yflux(i,j,k)     ...
                -0.25*(aam(i,j,k)+aam(i,j-1,k)     ...
                +aam(i,j,k-1)+aam(i,j-1,k-1))     ...
                *(h(i,j)+h(i,j-1))     ...
                *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)     ...
                /(dy(i,j)+dy(i,j-1));
            xflux(i,j,k)=0.5*(dy(i,j)+dy(i-1,j))*xflux(i,j,k);
            yflux(i,j,k)=0.5*(dx(i,j)+dx(i,j-1))*yflux(i,j,k);
        end
    end
end
%
%     do vertical advection: add flux terms, then step forward in time:
%
for k=2:kbm1
    for j=2:jmm1
        for i=2:imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))     ...
                *art(i,j)/(dz(k)+dz(k-1))     ...
                +xflux(i+1,j,k)-xflux(i,j,k)     ...
                +yflux(i,j+1,k)-yflux(i,j,k);
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)     ...
                *qb(i,j,k)-dti2*qf(i,j,k))     ...
                /((h(i,j)+etf(i,j))*art(i,j));
        end
    end
end
return 


