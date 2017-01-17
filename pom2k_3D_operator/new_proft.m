function [f1] = new_proft(f,wfsurf,fsurf,nbc,h,etf,swrad,kh)
%-----------------------------------------------------------------------
%     Irradiance parameters after Paulson and Simpson:
%       ntp               1      2       3       4       5
%   Jerlov type           i      ia      ib      ii     iii
r=[0.58,0.62,0.67,0.77,0.78];
ad1=[0.35,0.60,1.0,1.5,1.4];
ad2=[0.23,20.0,17.0,14.0,7.90];
% solves the equation:dti2*(kh*f')'-f=-fb
global im  jm kb dz_3d dzz_3d kbm1 dti2 umol kbm2 dz

rad=zeros(im,jm,kb);
a = zeros(im,jm,kb);c = zeros(im,jm,kb); d=zeros(im,jm,kb);
dh = h+etf;
dh_3d=repmat(dh,1,1,kb);swrad_3d=repmat(swrad,1,1,kb);
la=zeros(kb); f1=zeros(im,jm,kb);

if(nbc==2||nbc==4)
   rad=swrad_3d.*(r(ntp)*exp(z_3d.*dh_3d/ad1(ntp))+(1.e0-r(ntp))*exp(z_3d.*dh_3d/ad2(ntp)));
   rad(:,:,kb)=0.0;
end

    a(:,:,1:kbm2)=-dti2*(kh(:,:,2:kbm1)+umol); 
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);

    c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    c=DIVISION(-dti2*(kh+umol),dz_3d.*c.*dh_3d.*dh_3d);

    d=-f - DIVISION(dti2 .* DZF(rad),(dh_3d .* dz_3d));
    d(:,:,1)= -f(:,:,1) + dti2 .* wfsurf(:,:) ./ (dh(:,:) .* dz(1));
    d(:,:,kb)=-f(:,:,kb);

if(nbc==2)
    d(:,:,1)=-f(:,:,1) + dti2*(wfsurf+rad(:,:,1)-rad(:,:,2))./(dz(1)*dh);
elseif(nbc==3 || nbc==4)
    a(:,:,1)=0.e0;
    d(:,:,1)=-fsurf;
end
    
%     calculate penetrative radiation. At the bottom any unattenuated
%     radiation is deposited in the bottom layer:
    temp1=a(:,:,1:kbm1);    temp2=c(:,:,2:kb);  

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,:)+c(i,j,:)-1,kb,1),0) ...
      - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
      - diag(reshape(temp2(i,j,:),kbm1,1),-1);
   f1(i,j,:)=la\reshape(d(i,j,:),kb,1); 
      end
  end
return
