function [uf,wubot] =cal_u(advx,dt_3d,e_atmos,drhox,ub,u,v,w,egf,egb,etf,etb,km,vb,wusurf,cbc)
global im kb dti2 grav h_3d jm dz_3d dzz_3d kbm1 umol kbm2 dz dum_3d gs h cor;

uf=DIVISION( (AXB(repmat(etb,1,1,kb)+h_3d).*ub -dti2*( advx + drhox - AXB( repmat(cor,1,1,kb) .*dt_3d.*AYF(v) )...
   +grav*AXB(dt_3d).*( DXB( repmat(egf+egb,1,1,kb) )+DXB( repmat(e_atmos,1,1,kb) )*2.0 )/2.e0-DZF(AXB(w) .* AZB(u)))), ...
   AXB( repmat(etf,1,1,kb)+h_3d ) );
uf(:,:,kb)=0.e0;  bond=AXB(w) .* AZB(u); uf(im,:,:) = bond(im,:,:) ;   %add by hx

dh = AXB(h+etf);dh(1,:)=1.e0;dh(:,1)=1.e0;
dh_3d=repmat(dh,1,1,kb);
a = create_field(zeros(im,jm,kb),gs,6);
c = AXB(km);
d = create_field(zeros(im,jm,kb),gs,2);

    a(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);

    d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);
   
    tps = AXB(cbc) .* sqrt( ub(:,:,kbm1).^2 + AXB( AYF( vb(:,:,kbm1) ) ).^2 );
    a(:,:,kbm1)=-tps(:,:) * dti2./(dz(kbm1) .* dh(:,:));

    d=-uf;
    d(:,:,1)= -uf(:,:,1) + dti2 .* wusurf(:,:) ./ (dh(:,:) .* dz(1));
    temp1=a(:,:,1:kbm1);    temp1(:,:,kbm1)=0.e0;
    temp2=c(:,:,2:kb);

  for j=2:jm
      for i=2:im
        la=diag(reshape(a(i,j,:)+c(i,j,:)-1,kb,1),0) ...
	   - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
	   - diag(reshape(temp2(i,j,:),kbm1,1),-1);
	uf(i,j,:)=la\reshape(d(i,j,:),kb,1);  
      end
  end
   uf=uf.*dum_3d;
   wubot=-tps.*uf(:,:,kbm1);
   
   [uf] = bcond3_u(uf,u);
   
return

