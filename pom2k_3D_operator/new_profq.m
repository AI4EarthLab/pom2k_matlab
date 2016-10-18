function [sm,sh,dh,cc,ee,gg,l,kq,km,kh,...
         uf,vf,q2b,q2lb,a,c]=profq(sm,sh,dh,cc,ee,gg,l,kq,km,kh,uf,vf,q2,q2b,q2lb,a,c,h,etf,rho,u,v,dt,wusurf,wubot,wvsurf,wvbot,t,s)
load('grid.mat');load('depth.mat');load('para.mat');load('operator.mat');load('xyz.mat');load('masks.mat');     
a = zeros(im,jm,kb);c = zeros(im,jm,kb);
l0=zeros(im,jm);kn = zeros(im,jm,kb);boygr=zeros(im,jm,kb);
gh=zeros(im,jm,kb);stf=zeros(im,jm,kb);      
dh = h + etf;
dt_3d=repmat(dt,1,1,kb);dh_3d=repmat(dh,1,1,kb);h_3d=repmat(h,1,1,kb);
      
%%%%%%---------------------------------------------------      
temp = zeros(im,jm,kb);
dzz2_3d = zeros(im,jm,kb);
 for k=2:kbm1
     temp(:,:,k)=dzz(k-1)*dz(k-1); 
     dzz2_3d(:,:,k)=dzz(k-1);
 end

a = -dti2*(AZF2(kq)+umol) .* DIVISION(1.0,dzz2_3d.*dz_3d.*(dh_3d.^2));
c = -dti2*(AZB2(kq)+umol) .* DIVISION(1.0,temp.*(dh_3d.^2));
a(:,:,kb)=0;c(:,:,kb)=0;
%
%-----------------------------------------------------------------------
%
%     The following section solves the equation:
%
%       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
%
%     Surface and bottom boundary conditions:
%
 const1=(16.6e0^(2.e0/3.e0))*sef;
     
%
% initialize fields that are not calculated on all boundaries
% but are later used there
  ee(:,:,1) = OP_L_XY*ee(:,:,1)*OP_R_XY;
  gg(:,:,1) = OP_L_XY*gg(:,:,1)*OP_R_XY;
  l0= OP_L_XY*l0*OP_R_XY;

%
utau2=zeros(im,jm);
utau2 = sqrt(AXF1_XY(wusurf).^2+AYF1_XY(wvsurf).^2);

% Wave breaking energy- a variant of Craig & Banner (1994)
% see Mellor and Blumberg, 2003.
ee(:,:,1)=0.e0;
gg(:,:,1)=(15.8*cbcnst)^(2./3.)*utau2;
% Surface length scale following Stacey (1999).
l0 = surfl.*utau2/grav;
%
uf(:,:,kb) = sqrt(AXF1_XY(wubot).^2+AYF1_XY(wvbot).^2)*const1;
%
%    Calculate speed of sound squared:
tp = t+tbias;
sp = s+sbias;

%     Calculate pressure in units of decibars:
cc=zeros(im,jm,kb);
p=grav*rhoref*(-zz_3d.*h_3d)*1.e-4;
cc=1449.1e0+.00821e0*p+4.55e0*tp -.045e0*tp.^2+1.34e0*(sp-35.0e0);
cc=cc./sqrt((1.e0-.01642e0*p./cc).*(1.e0-0.40e0*p./cc.^2));
cc(:,:,kb)=0.e0;
%
%     Calculate buoyancy gradient:
%
%
q2b=abs(q2b);
q2lb=abs(q2lb);
boygr=-grav*DZB1(rho).*DIVISION( 1.0,(dzz2_3d.*h_3d) );
% *** NOTE: comment out next line if dens fores not include pressure     ... 
     +(grav^2)./AZB1(cc.^2); 
 

%

l=abs(q2lb./q2b);
l0_3d=repmat(l0,1,1,kb);
%[i,j,k]=ind2sub(size(z_3d),find(z_3d>-0.5));
%l([i, j, k]) = max(l([i, j, k]), kappa * l0_3d([i, j, k]));
idx = find(z_3d > -0.5);
l(idx) = max(l(idx), kappa * l0_3d(idx));
gh=l.^2.*boygr./q2b;
gh=min(gh,.028e0);

%
l(:,:,1)=kappa*l0;l(:,:,kb)=0.e0;
gh(:,:,1)=0.e0;gh(:,:,kb)=0.e0;

%
%    Calculate production of turbulent kinetic energy:
%
kn=km*sef.*(AXF2( DZB1(u) ).^2+AYF2( DZB1(v) ).^2).*DIVISION( 1.e0,(dzz2_3d.*dh_3d).^2 )...
%   Add shear due to internal wave field     ...
             -shiw*km.*boygr;
kn=kn+kh.*boygr;

%%%%%% boundary %%%%%%%%%%%%%
kn(1,:,:)=0.e0;kn(im,:,:)=0.e0;


%
%  NOTE: Richardson # dep. dissipation correction (Mellor: 2001; Ezer, 2000),
%  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
ghc=-6.0e0;
stf=ones(im,jm,kb);
      % It is unclear yet if diss. corr. is needed when surf. waves are included.
%           if(gh(i,j,k).lt.0.e0)
%    ...        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
%           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
dtef=sqrt(abs(q2b)).*stf./(b1*l+small);
%
%%%%% gg depends on ee during cycle
for k=2:kbm1
    gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-(2.e0*dti2*dtef(:,:,k)+1.e0));
    ee(:,:,k)=a(:,:,k).*gg(:,:,k);
    gg(:,:,k)=(-2.e0*dti2*kn(:,:,k)+c(:,:,k).*gg(:,:,k-1)-uf(:,:,k)).*gg(:,:,k);    
end

%
for k=1:kbm1
        ki=kb-k;
        uf(:,:,ki)=ee(:,:,ki).*uf(:,:,ki+1)+gg(:,:,ki);
end
%
%-----------------------------------------------------------------------
%
%     The following section solves the equation:
%
%       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
%
ee(:,:,2)=0.e0;gg(:,:,2)=0.e0;vf(:,:,kb)=0.e0;
dtef=dtef.*(1.e0+e2*((1.e0/abs(z_3d-z(1))+1.e0/abs(z_3d-z(kb))).*l1./(dh_3d*kappa)).^2);          
         
for k=2:kbm1
    gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-(dti2*dtef(:,:,k)+1.e0));
    ee(:,:,k)=a(:,:,k).*gg(:,:,k);
    gg(:,:,k)=(dti2*(-kn(:,:,k).*l1(:,:,k)*e1)+c(:,:,k).*gg(:,:,k-1)-vf(:,:,k)).*gg(:,:,k); 
end

%%%------------------------------------------------

for k=1:kb-2
    ki=kb-k;
    vf(:,:,ki)=ee(:,:,ki).*vf(:,:,ki+1)+gg(:,:,ki);
end
%
%%%%%%%%% boundary is not equal %%%%%%%%%%%%%%%%%%%
idx=find(uf<=small | vf<=small);
uf(idx)=small;
vf(idx)=0.1*dt_3d(idx)*small;
  
%
% % %-----------------------------------------------------------------------
%
%     The following section solves for km and kh:
%
      coef4=18.e0*a1*a1+9.e0*a1*a2;
      coef5=9.e0*a1*a2;
 %
%     Note that sm and sh limit to infinity when gh approaches 0.0288:
%
            coef1=a2*(1.e0-6.e0*a1/b1.*stf);
            coef2=3.e0*a2*b2./stf+18.e0*a1*a2;
            coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf);
            sh=coef1./(1.e0-coef2.*gh);
            sm=coef3+sh.*coef4.*gh;
            sm=sm./(1.e0-coef5.*gh);            
%
            kn=l.*sqrt(abs(q2));
            kq=(kn*.41e0.*sh+kq)*.5e0;
            km=(kn.*sm+km)*.5e0;
            kh=(kn.*sh+kh)*.5e0;
% %  
% cosmetics: make boundr. values as interior
% (even if not used: printout otherwise may show strange values)
           km(:,jm,:)=permute(km(:,jmm1,:),[1,3,2]) .* repmat(fsm(:,jm),1,kb);
           kh(:,jm,:)=permute(kh(:,jmm1,:),[1,3,2]) .* repmat(fsm(:,jm),1,kb);
           km(:,1,:)=permute(km(:,2,:),[1,3,2]) .* repmat(fsm(:,1),1,kb);
           kh(:,1,:)=permute(kh(:,2,:),[1,3,2]) .* repmat(fsm(:,1),1,kb);
           km(im,:,:)=permute(km(imm1,:,:),[2,3,1]) .* repmat(fsm(im,:)',1,kb);
           kh(im,:,:)=permute(kh(imm1,:,:),[2,3,1]) .* repmat(fsm(im,:)',1,kb);
           km(1,:,:)=permute(km(2,:,:),[2,3,1]) .* repmat(fsm(1,:)',1,kb);
           kh(1,:,:)=permute(kh(2,:,:),[2,3,1]) .* repmat(fsm(1,:)',1,kb);

return