function [uf] = advu(advx,aru,dz,cor,dt,e_atmos,dy,drhox,h,dti2,ub,uf,u,v,w,im,jm,kb,imm1,jmm1,kbm1,grav,egf,egb,etf,etb)
% **********************************************************************
% *                                                                    *
% * ROUTINE NAME:  advu                                                *
% *                                                                    *
% * FUN%TION    :  fores horizontal and vertical advection of           *
% *                u-momentum, and includes coriolis, surface slope    *
% *                and baroclinic terms.                               *
% *                                                                    *
% **********************************************************************
%
%
%     for vertical advection:
%


uf=zeros(im,jm,kb);
%
for k=2:kbm1
    for j=1:jm
        for i=2:im
            uf(i,j,k)=.25e0*(w(i,j,k)+w(i-1,j,k))  ...
                *(u(i,j,k)+u(i,j,k-1));
        end
    end
end
%
%     %ombine horizontal and vertical advection with coriolis, surface
%     slope and baroclinic terms:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            uf(i,j,k)=advx(i,j,k)     ...
                +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)     ...
                -aru(i,j)*.25e0     ...
                *(cor(i,j)*dt(i,j)     ...
                *(v(i,j+1,k)+v(i,j,k))     ...
                +cor(i-1,j)*dt(i-1,j)     ...
                *(v(i-1,j+1,k)+v(i-1,j,k)))     ...
                +grav*.125e0*(dt(i,j)+dt(i-1,j))     ...
                *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)     ...
                +(e_atmos(i,j)-e_atmos(i-1,j))*2.e0)     ...
                *(dy(i,j)+dy(i-1,j))     ...
                +drhox(i,j,k);
        end
    end
end
%
%     Step forward in time:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j)) ...
                *aru(i,j)*ub(i,j,k)     ...
                -2.e0*dti2*uf(i,j,k))     ...
                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))     ...
                *aru(i,j));
        end
    end
end

return
end
%


