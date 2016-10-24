function [uf,wubot] = new_profu(uf,etf,h,km,wusurf,cbc,ub,vb)

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
load('grid.mat');load('operator.mat');load('para.mat');
dh = AXB(h+etf);dh(1,:)=1.e0;dh(:,1)=1.e0;
a=zeros(im,jm,kb);ee=zeros(im,jm,kb);gg=zeros(im,jm,kb);
c = AXB(km);
%
for k=2:kbm1
    a(:,:,k-1) = DIVISION( -dti2*(c(:,:,k)+umol),(dz(k-1)*dzz(k-1)*dh.^2) );
    c(:,:,k) = DIVISION( -dti2*(c(:,:,k)+umol),(dz(k)*dzz(k-1)*dh.^2) );
end
%
ee(:,:,1)=a(:,:,1) ./ (a(:,:,1)-1.e0);
gg(:,:,1)=(-dti2*wusurf./(-dz(1)*dh)-uf(:,:,1)) ./ (a(:,:,1)-1.e0);
%
for k=2:kbm2
    gg(:,:,k)=1.e0 ./ (a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
    ee(:,:,k)=a(:,:,k).* gg(:,:,k);
    gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-uf(:,:,k)).*gg(:,:,k);
end
%
tps = AXB(cbc) .* sqrt( ub(:,:,kbm1).^2 + AXB( AYF( vb(:,:,kbm1) ) ).^2 );
uf(:,:,kbm1) = (c(:,:,kbm1).* gg(:,:,kbm2)-uf(:,:,kbm1))...
               ./(tps*dti2./(-dz(kbm1)*dh)-1.e0-(ee(:,:,kbm2)-1.e0).*c(:,:,kbm1)).*dum;
%
for k=2:kbm1
    ki=kb-k;
    uf(:,:,ki)=(ee(:,:,ki).*uf(:,:,ki+1)+gg(:,:,ki)).*dum;
end
%
wubot=-tps.*uf(:,:,kbm1);

return

