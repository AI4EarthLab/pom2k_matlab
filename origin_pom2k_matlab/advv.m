function [vf] = advv(advy,arv,dz,cor,dt,e_atmos,dx,drhoy,h,dti2,vb,vf,u,v,w,im,jm,kb,imm1,jmm1,kbm1,grav,egf,egb,etf,etb)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  fores horizontal and vertical advection of           *
% *                v-momentum, and includes coriolis, surface slope    *
% *                and baroclinic terms.                               *
% *                                                                    *
% **********************************************************************
%
%
%     for vertical advection:
%


vf=zeros(im,jm,kb);

for k=2:kbm1
    for j=2:jm
        for i=1:im
            vf(i,j,k)=.25e0*(w(i,j,k)+w(i,j-1,k))     ...
                *(v(i,j,k)+v(i,j,k-1));
        end
    end
end
%
%     %ombine horizontal and vertical advection with coriolis: surface
%     slope and baroclinic terms:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            vf(i,j,k)=advy(i,j,k)     ...
                +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)     ...
                +arv(i,j)*.25e0     ...
                *(cor(i,j)*dt(i,j)     ...
                *(u(i+1,j,k)+u(i,j,k))     ...
                +cor(i,j-1)*dt(i,j-1)     ...
                *(u(i+1,j-1,k)+u(i,j-1,k)))     ...
                +grav*.125e0*(dt(i,j)+dt(i,j-1))     ...
                *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)     ...
                +(e_atmos(i,j)-e_atmos(i,j-1))*2.e0)     ...
                *(dx(i,j)+dx(i,j-1))     ...
                +drhoy(i,j,k);
        end
    end
end
%
%     Step forward in time:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))     ...
                *arv(i,j)*vb(i,j,k)     ...
                -2.e0*dti2*vf(i,j,k))     ...
                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))     ...
                *arv(i,j));
        end
    end
end
%
return
%

end
%

