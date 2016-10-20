function [f] = new_proft(f,wfsurf,fsurf,nbc,h,etf,swrad,kh)

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
load('grid.mat');load('operator.mat');load('para.mat');

dh = h+etf;
rad=zeros(im,jm,kb);
%
a = zeros(im,jm,kb);c = zeros(im,jm,kb);
ee = zeros(im,jm,kb);gg = zeros(im,jm,kb);
dh = h + etf;
dh_3d=repmat(dh,1,1,kb);swrad_3d=repmat(swrad,1,1,kb);

for k=2:kbm1
    a(:,:,k-1)=-dti2*(kh(:,:,k)+umol)./(dz(k-1)*dzz(k-1).*dh(:,:).*dh(:,:));
    c(:,:,k)  =-dti2*(kh(:,:,k)+umol)./(dz(k)  *dzz(k-1).*dh(:,:).*dh(:,:));
end

%
%     calculate penetrative radiation. At the bottom any unattenuated
%     radiation is deposited in the bottom layer:
if(nbc==2||nbc==4)
   rad=swrad_3d.*(r(ntp)*exp(z_3d.*dh_3d/ad1(ntp))+(1.e0-r(ntp))*exp(z_3d.*dh_3d/ad2(ntp)));
   rad(:,:,kb)=0.0;
end

if(nbc==1)
    ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
    gg(:,:,1)=DIVISION( (-dti2*wfsurf./(-dz(1)*dh)-f(:,:,1)),a(:,:,1)-1.e0 );
elseif(nbc==2)
    ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
    gg(:,:,1)=DIVISION( (dti2*(wfsurf+rad(:,:,1)-rad(:,:,2))./(dz(1)*dh)-f(:,:,1)),a(:,:,1)-1.e0 );
elseif(nbc==3 || nbc==4)
    ee(:,:,1) = 0.e0;
    gg(:,:,1) = fsurf;
end

for k=2:kbm2
    gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
    ee(:,:,k)=a(:,:,k).*gg(:,:,k);
    gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-f(:,:,k)+dti2*(rad(:,:,k)-rad(:,:,k+1))./(dh*dz(k))).*gg(:,:,k);
end

% Bottom adiabatic boundary condition:
f(:,:,kbm1)=(c(:,:,kbm1).*gg(:,:,kbm2)-f(:,:,kbm1) + dti2*(rad(:,:,kbm1)-rad(:,:,kb))./(dh*dz(kbm1))) ...
            ./(c(:,:,kbm1).*(1.e0-ee(:,:,kbm2))-1.e0);
        
for k=2:kbm1
    ki=kb-k;
    f(:,:,ki)=(ee(:,:,ki).*f(:,:,ki+1)+gg(:,:,ki));
end

return
