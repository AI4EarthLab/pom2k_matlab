function [a,c,ee,gg,tps,uf,wubot] = profu(a,c,ee,gg,tps,uf,wubot,...
                                    etf,h,km,dti2,umol,dz,dzz,wusurf,cbc,dum,im,jm,kb,imm1,jmm1,kbm1,kbm2,ub,vb)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  Solves for vertical diffusion of x-momentum using   *
% *                method described by Richmeyer and Morton.           *
% *                                                                    *
% *                See:                                                *
% *                                                                    *
% *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
% *                  Methods for Initial-Value Problems, 2nd edition,  *
% *                  Interscience, New York, 198-201.                  *
% *                                                                    *
% *                NOTE that wusurf has the opposite sign to the wind  *
% *                speed.                                              *
% *                                                                    *
% **********************************************************************
%
dh=ones(im,jm);
%
%     The following section solves the equation:
%
%       dti2*(km*u')'-u=-ub
%
%
for j=2:jm
    for i=2:im
        dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5e0;
    end
end
%
for k=1:kb
    for j=2:jm
        for i=2:im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5e0;
        end
    end
end
%
for k=2:kbm1
    for j=1:jm
        for i=1:im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol) ...
                /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j));
            c(i,j,k)=-dti2*(c(i,j,k)+umol)     ...
                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j));
        end
    end
end
%
for j=1:jm
    for i=1:im
        ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0);
        gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))     ...
            -uf(i,j,1))     ...
            /(a(i,j,1)-1.e0);
    end
end
%
for k=2:kbm2
    for j=1:jm
        for i=1:im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0);
            ee(i,j,k)=a(i,j,k)*gg(i,j,k);
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k);
        end
    end
end
%
for j=2:jmm1
    for i=2:imm1
        tps(i,j)=0.5e0*(cbc(i,j)+cbc(i-1,j))     ...
            *sqrt(ub(i,j,kbm1)^2     ...
            +(.25e0*(vb(i,j,kbm1)+vb(i,j+1,kbm1)     ...
            +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))^2);
        uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))     ...
            /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0     ...
            -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1));
        uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j);
    end
end
%

for k=2:kbm1
    ki=kb-k;
    for j=2:jmm1
        for i=2:imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j);
        end
    end
end
%
for j=2:jmm1
    for i=2:imm1
        wubot(i,j)=-tps(i,j)*uf(i,j,kbm1);
    end
end

return 
end

%
