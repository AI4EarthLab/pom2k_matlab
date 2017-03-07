 function [l,kq,km,kh,uf,vf,q2b,q2lb]=profq(kq,km,kh,uf,...
           vf,q2,q2b,q2lb,etf,rho,u,v,dt,wusurf,wubot,wvsurf,wvbot,t,s)

global im jm kb imm1 jmm1 kbm1 kbm2 fsm_3d dzz_3d dz_3d zz_3d z_3d gs ...
       small umol dti2 grav kappa tbias sbias rhoref z h h_3d

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
%-----------------------------------------------------------------------
%     The following section solves the equation:
%
%       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
%
%     Surface and bottom boundary conditions:
%
const1=(16.6e0^(2.e0/3.e0))*sef;      
utau2 = sqrt( AXF(wusurf).^2 +AYF(wvsurf).^2 );
l0 = surfl*utau2/grav;                     
uf(:,:,kb) = sqrt( AXF(wubot).^2 +AYF(wvbot).^2 ) .* const1;
        
%    Calculate speed of sound squared:    
p=grav*rhoref*(-zz_3d .* h_3d)*1.e-4;     %     Calculate pressure in units of decibars:
cc=1449.10+.00821*p+4.55*(t+tbias) -.045e0*(t+tbias).^2 +1.34*(s+sbias-35.0e0);
cc=cc./sqrt((1.0-.01642.*p./cc) .*(1.0-0.40.*p./cc.^2));      
cc(:,:,kb)=0;
      
%     Calculate buoyancy gradient:
q2b =abs(q2b);
q2lb=abs(q2lb);
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
%
%  NOTE: Richardson # dep. dissipation correction (Mellor: 2001; Ezer, 2000),
%  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
ghc=-6.0e0;
stf=ones(im,jm,kb);
dtef=sqrt(q2b).*stf./(b1.*l+small);
dtef(:,:,1)=0.e0;dtef(:,:,kb)=0.e0;

    d=-uf-2.e0*dti2*kn;
    d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
    d(:,jm,1)=0; d(im,:,1)=0; d(:,:,kb)=-uf(:,:,kb);
    temp1=a(:,:,1:kbm1);    temp2=c(:,:,2:kb);

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,:)+c(i,j,:)-1-2.e0*dti2.*dtef(i,j,:),kb,1),0) ...
      - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
      - diag(reshape(temp2(i,j,:),kbm1,1),-1);
   uf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
      end
  end
    
%-----------------------------------------------------------------------
%     The following section solves the equation:
%
%       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb

vf(:,:,kb)=0.e0;

for k=2:kbm1
    dtef(:,:,k)=dtef(:,:,k).*(1.e0+e2.*((1.e0/abs(z(k)-z(1))+1.e0./abs(z(k)-z(kb)))     ...
                        .*l(:,:,k)./(dh(:,:)*kappa)).^2);
end

    d=-vf-dti2*kn.*l.*e1;
    d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
    d(:,jm,1)=0.e0; d(im,:,1)=0.e0; d(:,:,kb)=0.e0;
    temp=vf(:,:,1);

  for j=2:jm
      for i=2:im
   la=diag(reshape(a(i,j,:)+c(i,j,:)-1-dti2.*dtef(i,j,:),kb,1),0) ...
      - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
      - diag(reshape(temp2(i,j,:),kbm1,1),-1);
   vf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
      end
  end
   vf(:,:,1)=temp;

dt_3d=repmat(dt,1,1,kb);

filter = (uf<=small | vf<=small);
filter(:,:,1) = false;
filter(:,:,kb) = false;
uf(filter.data) = small;
vf(filter.data) = 0.1 * dt_3d(filter.data) * small;

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
      
 end
 
