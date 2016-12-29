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
%load('grid.mat');load('operator.mat');load('para.mat');
global im  jm kb dz_3d dzz_3d kbm1 dti2 umol kbm2 dz dum

dh = AXB(h+etf);dh(1,:)=1.e0;dh(:,1)=1.e0;
dh_3d=repmat(dh,1,1,kb);
a=zeros(im,jm,kb);ee=zeros(im,jm,kb);gg=zeros(im,jm,kb);
c = AXB(km); c(1,:,:)=0.e0; c(:,1,:)=0.e0;
la=zeros(kbm1);d=zeros(im,jm,kb);

    d(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
    a=DIVISION(d,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
    a(:,:,kbm1)= 0.e0;

    d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);
    c(:,:,1)  =0.e0;
   
tps = AXB(cbc) .* sqrt( ub(:,:,kbm1).^2 + AXB( AYF( vb(:,:,kbm1) ) ).^2 );

    d=-uf;
    d(:,:,1)= -uf(:,:,1) + dti2 .* wusurf(:,:) ./ (dh(:,:) .* dz(1));
    d(:,:,kbm1)=-uf(:,:,kbm1)+tps.*dti2./(dz(kbm1)*dh); 

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,1:kbm1)+c(i,j,1:kbm1)-1,kbm1,1),0) ...
      - diag(reshape(a(i,j,1:kbm2),kbm2,1),1) ...
      - diag(reshape(c(i,j,2:kbm1),kbm2,1),-1);
   uf(i,j,1:kbm1)=la\reshape(d(i,j,1:kbm1),kbm1,1); 
      end
  end
   uf=uf.*repmat(dum,1,1,kb);
   wubot=-tps.*uf(:,:,kbm1);

return

