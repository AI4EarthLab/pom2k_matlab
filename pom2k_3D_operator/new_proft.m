function [f,wfsurf,fsurf,nbc,dh,...
         a,c,ee,gg] = proft(f,wfsurf,fsurf,nbc,dh,...
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
load('grid.mat');load('operator.mat');load('para.mat');
a = zeros(im,jm,kb);c = zeros(im,jm,kb);
ee = zeros(im,jm,kb);gg = zeros(im,jm,kb);
tmp1 = zeros(im,jm,kb);tmp2 = zeros(im,jm,kb);rad=zeros(im,jm,kb);
dh = h + etf;
dh_3d=repmat(dh,1,1,kb);swrad_3d=repmat(swrad,1,1,kb);

for k=2:kbm1
    tmp1(:,:,k)=dzz(k-1)*dz(k);
    tmp2(:,:,k)=dzz(k-1)*dz(k-1);
end
aa= DIVISION(-dti2 .* (kh+umol) , (tmp2 .* dh_3d.^2)); 
c= DIVISION(-dti2 .* (kh+umol) , (tmp1 .* dh_3d.^2));
a(:,:,1:kbm2)=aa(:,:,2:kbm1);

if(nbc==2||nbc==4)
   rad=swrad_3d.*(r(ntp)*exp(z_3d.*dh_3d/ad1(ntp))...
       +(1.e0-r(ntp))*exp(z_3d.*dh_3d/ad2(ntp)));
   rad(:,:,kb)=0.e0;
end
    %
if(nbc==1)
   ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
   gg(:,:,1)=DIVISION( (-dti2*wfsurf./(-dz(1)*dh)-f(:,:,1)),a(:,:,1)-1.e0 );
    %
elseif(nbc==2)
    ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
    gg(:,:,1)=DIVISION( (dti2*(wfsurf+rad(:,:,1)-rad(:,:,2))./(dz(1)*dh)-f(:,:,1)),a(:,:,1)-1.e0 );
    %
elseif(nbc==3 || nbc==4)
    ee(:,:,1) = 0.e0;
    gg(:,:,1) = fsurf;
end

for k=2:kbm2
    gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
    ee(:,:,k)=a(:,:,k).*gg(:,:,k);
    gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-f(:,:,k)+dti2*(rad(:,:,k)-rad(:,:,k+1))...
                ./(dh*dz(k))).*gg(:,:,k);
end
%
%     Bottom adiabatic boundary condition:
%
f(:,:,kbm1)=(c(:,:,kbm1).*gg(:,:,kbm2)-f(:,:,kbm1) ...
            +dti2*(rad(:,:,kbm1)-rad(:,:,kb))./(dh*dz(kbm1))) ...
            ./(c(:,:,kbm1).*(1.e0-ee(:,:,kbm2))-1.e0);

for k=2:kbm1
    ki=kb-k;
    f(:,:,ki)=(ee(:,:,ki).*f(:,:,ki+1)+gg(:,:,ki));
end
return
