function [f,wfsurf,fsurf,nbc,dh,a,c,ee,gg] = proft(f,wfsurf,fsurf,nbc,dh,...
                            a,c,ee,gg,...
                            h,etf,dti2,dz,dzz,swrad,ntp,im,jm,kb,kbm1,kbm2,kh,umol)

% **********************************************************************
% *                                                                    *
% * FUN%TION    :  Solves for vertical diffusion of temperature and    *
% *                salinity using method described by Richmeyer and    *
% *                Morton.                                             *
% *                                                                    *
% *                Irradiance parameters are from Paulson and Simpson. *
% *                                                                    *
% *                See:                                                *
% *                                                                    *
% *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
% *                  Methods for Initial-Value Problems, 2nd edition,  *
% *                  Interscience, New York, 198-201.                  *
% *                                                                    *
% *                Paulson, %. A., and J. Simpson, 1977: Irradiance    *
% *                  measurements in the upper ocean, J. Phys.         *
% *                  Oceanogr., 7, 952-956.                            *
% *                                                                    *
% *                NOTES:                                              *
% *                                                                    *
% *                (1) wfsurf and swrad are negative values when water *
% *                    column is warming or salt is being added.       *
% *                                                                    *
% *                (2) nbc may only be 1 and 3 for salinity.           *
% *                                                                    *
% **********************************************************************
%
%

%
%-----------------------------------------------------------------------
%
%     Irradiance parameters after Paulson and Simpson:
%
%       ntp               1      2       3       4       5
%   Jerlov type           i      ia      ib      ii     iii
%
r=[0.58,0.62,0.67,0.77,0.78];
ad1=[0.35,0.60,1.0,1.5,1.4];
ad2=[0.23,20.0,17.0,14.0,7.90];
%
%-----------------------------------------------------------------------
%
%     Surface boundary condition:
%
%       nbc   prescribed    prescribed   short wave
%             temperature      flux      penetration
%             or salinity               (temperature
%                                           only)
%
%        1        no           yes           no
%        2        no           yes           yes
%        3        yes          no            no
%        4        yes          no            yes
%
%     NOTE that only 1 and 3 are allowed for salinity.
%
%-----------------------------------------------------------------------
%
%     The following section solves the equation:
%
%       dti2*(kh*f')'-f=-fb
%

dh = h+etf;
rad=zeros(im,jm,kb);
%

for k=2:kbm1
    for j=1:jm
        for i=1:im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)  ...
                /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j));
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)     ...
                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j));
        end
    end
end
%
%     %alculate penetrative radiation. At the bottom any unattenuated
%     radiation is deposited in the bottom layer:
%
%
if(nbc==2||nbc==4)
    %
    for k=1:kbm1
        for j=1:jm
            for i=1:im
                rad(i,j,k)=swrad(i,j)     ...
                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp))     ...
                    +(1.e0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)));
            end
        end
    end
    %
end
%
if(nbc==1)
    %
    for j=1:jm
        for i=1:im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0);
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1);
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0);
        end
    end
    %
elseif(nbc==2)
    %
    for j=1:jm
        for i=1:im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0);
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))     ...
                /(dz(1)*dh(i,j))     ...
                -f(i,j,1);
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0);
        end
    end
    %
elseif(nbc==3 || nbc==4)
    %
    for j=1:jm
        for i=1:im
            ee(i,j,1)=0.e0;
            gg(i,j,1)=fsurf(i,j);
        end
    end
    %
end
%
for k=2:kbm2
    for j=1:jm
        for i=1:im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0);
            ee(i,j,k)=a(i,j,k)*gg(i,j,k);
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)     ...
                +dti2*(rad(i,j,k)-rad(i,j,k+1))     ...
                /(dh(i,j)*dz(k)))     ...
                *gg(i,j,k);
        end
    end
end
%
%     Bottom adiabatic boundary condition:
%
for j=1:jm
    for i=1:im
        f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)     ...
            +dti2*(rad(i,j,kbm1)-rad(i,j,kb))     ...
            /(dh(i,j)*dz(kbm1)))     ...
            /(c(i,j,kbm1)*(1.e0-ee(i,j,kbm2))-1.e0);
    end
end
%
for k=2:kbm1
    ki=kb-k;
    for j=1:jm
        for i=1:im
            f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki));
        end
    end
end

%
