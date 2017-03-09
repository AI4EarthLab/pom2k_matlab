function [tf]=cal_ts(tb,t,dt_3d,u,v,aam,w,etb,etf,wtsurf,tsurf,nbct,swrad,kh)
global kb dum_3d dvm_3d dti2 tprni h_3d im jm dz_3d dzz_3d kbm1 umol kbm2 dz gs h

tf=((h_3d+repmat(etb,1,1,kb)).*tb - dti2 .* (DXF( AXB(dt_3d).*AXB(t).*u-AXB(aam).*AXB(h_3d)*tprni.*DXB(tb).*dum_3d ) ...
   + DYF( AYB(dt_3d).*AYB(t).*v-AYB(aam).*AYB(h_3d)*tprni.*DYB(tb).*dvm_3d )-DZF( AZB(t).*w ))) ./((h_3d+repmat(etf,1,1,kb))) ;

%------------------------------------------------------------------
r=[0.58,0.62,0.67,0.77,0.78];
ad1=[0.35,0.60,1.0,1.5,1.4];
ad2=[0.23,20.0,17.0,14.0,7.90];

rad=create_field(zeros(im,jm,kb),gs,7);
a = create_field(zeros(im,jm,kb),gs,7);
c = create_field(zeros(im,jm,kb),gs,7);

dh = h+etf;
dh_3d=repmat(dh,1,1,kb);

if(nbct==2||nbct==4)
   rad=repmat(swrad,1,1,kb) .*(r(ntp)*exp(z_3d.*dh_3d/ad1(ntp))+(1.e0-r(ntp))*exp(z_3d.*dh_3d/ad2(ntp)));
   rad(:,:,kb)=0.0;
end

    a(:,:,1:kbm2)=-dti2*(kh(:,:,2:kbm1)+umol); 
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);

    c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    c=DIVISION(-dti2*(kh+umol),dz_3d.*c.*dh_3d.*dh_3d);

    d=-tf - dti2 .* DZF(rad)./dh_3d;
    d(:,:,1)= -tf(:,:,1) + dti2 .* wtsurf(:,:) ./ (dh(:,:) .* dz(1));

if(nbct==2)
    d(:,:,1)=-tf(:,:,1) + dti2*(wtsurf+rad(:,:,1)-rad(:,:,2))./(dz(1)*dh);
elseif(nbct==3 || nbct==4)
    a(:,:,1)=0.e0;     d(:,:,1)=-tsurf;
end
    
    temp1=a(:,:,1:kbm1);    temp2=c(:,:,2:kb);  

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,:)+c(i,j,:)-1,kb,1),0) ...
      - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
      - diag(reshape(temp2(i,j,:),kbm1,1),-1);
   tf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
      end
  end
return