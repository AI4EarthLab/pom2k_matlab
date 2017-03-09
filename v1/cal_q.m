function [q2f,q2,q2b,q2lf,q2l,q2lb,km,kq,kh]=cal_q(q2b,q2,q2lb,q2l,dt_3d,u,v,w,aam,etb,etf,rho, ...
                                            wusurf,wubot,wvsurf,wvbot,km,kq,kh,t,s)

global im jm kb imm1 jmm1 kbm1 kbm2 dum_3d dvm_3d dti2 h_3d dzz_3d dz_3d zz_3d z_3d ...
       gs small umol grav kappa tbias sbias rhoref h z fsm_3d smoth;

q2f= ( (h_3d+repmat(etb,1,1,kb)) .* q2b -dti2*( -DZB(AZF(w.*q2)) + DXF(AXB(q2) .* AXB(dt_3d) .* AZB(u) ... 
     -AZB( AXB( aam )).*AXB(h_3d) .*DXB( q2b ).* dum_3d) + DYF(AYB(q2) .* AYB(dt_3d) .* AZB(v) ...
     -AZB( AYB( aam )).*AYB(h_3d) .*DYB( q2b ).* dvm_3d)) )./  ( h_3d+repmat(etf,1,1,kb) ) ;

q2lf= ( (h_3d+repmat(etb,1,1,kb)) .* q2lb -dti2*( -DZB(AZF(w.*q2l)) + DXF(AXB(q2l) .* AXB(dt_3d) .* AZB(u) ...
    -AZB( AXB( aam )).*AXB(h_3d) .*DXB( q2lb ).* dum_3d) +DYF(AYB(q2l) .* AYB(dt_3d) .* AZB(v)  ...
    -AZB( AYB( aam )).*AYB(h_3d) .*DYB( q2lb ).* dvm_3d)) )./  ( h_3d+repmat(etf,1,1,kb) ) ;  
    
a1=0.92; b1=16.6 ; a2=0.74 ; b2=10.1    ; c1=0.08;
e1=1.8 ; e2=1.33 ; sef=1.0 ; cbcnst=100.; surfl=2.e5 ; shiw=0.0;
dh=h+etf;      
dh_3d=repmat(dh,1,1,kb);
a = create_field(zeros(im,jm,kb),gs,7);
c = create_field(zeros(im,jm,kb),gs,7);

a(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,2:kbm1);
c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,1:kbm2);

a= DIVISION(-dti2 .* (AZF(kq)+umol) , (a .* dh_3d .* dh_3d));
c= DIVISION(-dti2 .* (AZB(kq)+umol) , (c .* dh_3d .* dh_3d));

const1=(16.6e0^(2.e0/3.e0))*sef;      
utau2 = sqrt( AXF(wusurf).^2 +AYF(wvsurf).^2 );
l0 = surfl*utau2/grav;                     
q2f(:,:,kb) = sqrt( AXF(wubot).^2 +AYF(wvbot).^2 ) .* const1;
        
%    Calculate speed of sound squared:    
p=grav*rhoref*(-zz_3d .* h_3d)*1.e-4;     %     Calculate pressure in units of decibars:
cc=1449.10+.00821*p+4.55*(t+tbias) -.045e0*(t+tbias).^2 +1.34*(s+sbias-35.0e0);
cc=cc./sqrt((1.0-.01642.*p./cc) .*(1.0-0.40.*p./cc.^2));      
cc(:,:,kb)=0;
      
q2b =abs(q2b);  q2lb=abs(q2lb);
boygr=-grav* DZB(rho)./h_3d + DIVISION(grav^2 , AZB(cc.^2));
boygr(:,:,1)=0.e0;
l=q2lb ./ q2b;
%l=max(l, repmat(kappa*l0,1,1,kb));
tmp1=create_field((z_3d>-0.5),gs,7);
tmp2=create_field((z_3d<=-0.5),gs,7);
l=tmp1 .* max(l, repmat(kappa*l0,1,1,kb))+tmp2 .* l;
gh=l.^2 .* boygr ./q2b;
gh=min(gh, 0.028);

l(:,:,1)=kappa*l0; l(:,:,kb)=0;
gh(:,:,1)=0      ; gh(:,:,kb)=0;

%    Calculate production of turbulent kinetic energy:
kn= km.*sef.*(DZB(AXF(u)).^2 + DZB(AYF(v)).^2)./(dh_3d.^2) -shiw.*km.*boygr + kh.*boygr;

ghc=-6.0e0;
stf=ones(im,jm,kb);
dtef=sqrt(q2b).*stf./(b1.*l+small);
dtef(:,:,1)=0.e0;dtef(:,:,kb)=0.e0;

    d=-q2f-2.e0*dti2*kn;
    d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
    d(:,jm,1)=0; d(im,:,1)=0; d(:,:,kb)=-q2f(:,:,kb);
    temp1=a(:,:,1:kbm1);    temp2=c(:,:,2:kb);

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,:)+c(i,j,:)-1-2.e0*dti2.*dtef(i,j,:),kb,1),0) ...
      - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
      - diag(reshape(temp2(i,j,:),kbm1,1),-1);
   q2f(i,j,:)=la\reshape(d(i,j,:),kb,1); 
      end
  end
  
  %-----------------------------------------------------------------------
q2lf(:,:,kb)=0.e0;
for k=2:kbm1
    dtef(:,:,k)=dtef(:,:,k).*(1.e0+e2.*((1.e0/abs(z(k)-z(1))+1.e0./abs(z(k)-z(kb)))     ...
                        .*l(:,:,k)./(dh(:,:)*kappa)).^2);
end

    d=-q2lf-dti2*kn.*l.*e1;
    d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
    d(:,jm,1)=0.e0; d(im,:,1)=0.e0; d(:,:,kb)=0.e0;
    temp=q2lf(:,:,1);

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,:)+c(i,j,:)-1-dti2.*dtef(i,j,:),kb,1),0) ...
      - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
      - diag(reshape(temp2(i,j,:),kbm1,1),-1);
   q2lf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
      end
  end
   q2lf(:,:,1)=temp;

% dt_3d=repmat(dt,1,1,kb);

filter = (q2f<=small | q2lf<=small);
filter(:,:,1) = false;
filter(:,:,kb) = false;
q2f(filter.data) = small;
q2lf(filter.data) = 0.1 * dt_3d(filter.data) * small;

%-----------------------------------------------------------------------
%     The following section solves for km and kh:
    coef4=18.e0*a1*a1+9.e0*a1*a2;
    coef5=9.e0*a1*a2;
%     Note that sm and sh limit to infinity when gh approaches 0.0288:
    coef1=a2*(1.e0-6.e0*a1/b1*stf);
    coef2=3.e0*a2*b2/stf+18.e0*a1*a2;
    coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf);
    sh=coef1./(1.e0-coef2.*gh);
    sm=coef3+sh.*coef4.*gh;
    sm=sm./(1.0-coef5.*gh);
            
    kn=l.*sqrt(abs(q2));
    kq=(kn.*.41e0.*sh+kq)*.5e0;
    km=(kn.*sm+km)*.5e0;
    kh=(kn.*sh+kh)*.5e0;
       
      km(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
      kh(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
      km(:,1,:)=km(:,2,:).*fsm_3d(:,1,:);
      kh(:,1,:)=kh(:,2,:).*fsm_3d(:,1,:);
      
      km(im,:,:)=km(imm1,:,:).*fsm_3d(im,:,:);
      kh(im,:,:)=kh(imm1,:,:).*fsm_3d(im,:,:);
      km(1,:,:)=km(2,:,:).*fsm_3d(1,:,:);
      kh(1,:,:)=kh(2,:,:).*fsm_3d(1,:,:);

    [q2f, q2lf] = bcond6(q2f, q2lf, u, v, q2, q2l);
    [q2,q2l,q2b,q2lb]=smoth_update(q2f,q2lf,q2,q2l,q2b,q2lb);
 end