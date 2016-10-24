function [vf,wvbot] = new_profv(vf,etf,h,km,wvsurf,cbc,ub,vb)
% **********************************************************************
%                                                                      *
% * FUNCTION    :  Solves for vertical diffusion of y-momentum using   *
% *                method described by Richmeyer and Morton.           *
% *                                                                    *
% *                See:                                                *
% *                                                                    *
% *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
% *                  Methods for Initial-Value Problems, 2nd edition,  *
% *                  Interscience, New York, 198-201.                  *
% *                                                                    *
% *                NOTE that wvsurf has the opposite sign to the wind  *
% *                speed.                                              *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('operator.mat'); load('para.mat');
a=zeros(im,jm,kb);
ee=zeros(im,jm,kb); 
gg=zeros(im,jm,kb);
dh = AYB1(h+etf); dh(1,:)=1.e0; dh(:,1)=1.e0;
c = AYB1(km); c(1,:,:)=0.e0;

%
for k=2:kbm1
    a(:,:,k-1) = -dti2 * (c(:,:,k)+umol) ./ (dz(k-1)*dzz(k-1)*dh.^2);
    c(:,:,k) = -dti2 * (c(:,:,k)+umol) ./ (dz(k)*dzz(k-1)*dh.^2);
end
%
ee(:,:,1) = a(:,:,1) ./ (a(:,:,1)-1.e0);
gg(:,:,1) = (-dti2*wvsurf./(-dz(1)*dh)-vf(:,:,1)) ./ (a(:,:,1)-1.e0);
%
for k=2:kbm2
    gg(:,:,k) = 1.e0 ./ (a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
    ee(:,:,k) = a(:,:,k) .* gg(:,:,k);
    gg(:,:,k) = (c(:,:,k) .* gg(:,:,k-1)-vf(:,:,k)) .* gg(:,:,k);
end
%
tps = AYB2(cbc) .* sqrt( AYB2( AXF1( ub(:,:,kbm1) ) ).^2 + vb(:,:,kbm1).^2 );
tps(1,:) = 0.e0;
vf(:,:,kbm1) = (c(:,:,kbm1) .* gg(:,:,kbm2) - vf(:,:,kbm1))  ...
               ./(tps * dti2./(-dz(kbm1) * dh)-1.e0 - (ee(:,:,kbm2)-1.e0) .* c(:,:,kbm1)) .* dvm;
%
for k=2:kbm1
    ki=kb-k;
    vf(:,:,ki) = (ee(:,:,ki) .* vf(:,:,ki+1) + gg(:,:,ki)) .* dvm;
end
%
wvbot=-tps .* vf(:,:,kbm1);

return
