function [vf,wvbot] = cal_v(advy,dt_3d,e_atmos,drhoy,vb,u,v,w,egf,egb,etf,etb,km,ub,wvsurf,cbc)
global kb dti2 grav h_3d im jm dz_3d dzz_3d kbm1 umol kbm2 dz dvm_3d gs h cor;

vf = ( AYB( repmat(etb,1,1,kb) +h_3d).*vb-dti2*( advy+drhoy+AYB( repmat(cor,1,1,kb) .*dt_3d.*AXF(u) )...
     +grav*AYB(dt_3d).*( DYB( repmat(egf+egb,1,1,kb) )+DYB( repmat(e_atmos,1,1,kb) )*2.0 )/2.e0-DZF(AYB(w) .* AZB(v)))) ...
     ./AYB( repmat(etf,1,1,kb) +h_3d);
vf(:,:,kb)=0.e0;

dh = AYB(h+etf); dh(1,:)=1.e0; dh(:,1)=1.e0;
dh_3d=repmat(dh,1,1,kb);
a = create_field(zeros(im,jm,kb),gs,1);
c = AYB(km); 
d = create_field(zeros(im,jm,kb),gs,1);

    a(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);

    d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);

    tps = AYB(cbc) .* sqrt( AYB( AXF( ub(:,:,kbm1) ) ).^2 + vb(:,:,kbm1).^2 );
    a(:,:,kbm1)=-tps(:,:) * dti2./(dz(kbm1) .* dh(:,:));
		    
    d=-vf;
    d(:,:,1)= -vf(:,:,1) + dti2 .* wvsurf(:,:) ./ (dh(:,:) .* dz(1));    
    temp1=a(:,:,1:kbm1); temp1(:,:,kbm1)=0.e0; temp2=c(:,:,2:kb);
			        
    for j=2:jm
       for i=2:im
      la=diag(reshape(a(i,j,:)+c(i,j,:)-1,kb,1),0) ...
         - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
         - diag(reshape(temp2(i,j,:),kbm1,1),-1);
      vf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
        end
    end

   vf=vf.*dvm_3d;
   wvbot=-tps .* vf(:,:,kbm1);
return