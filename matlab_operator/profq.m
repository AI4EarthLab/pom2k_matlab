 function [sm,sh,dh,cc,ee,gg,l,kq,km,kh,...
         uf,vf,q2b,q2lb,a,c]=profq(sm,sh,dh,cc,...
         ee,gg,l,kq,km,kh,uf,vf,q2,q2b,q2lb,a,c,...
         h,etf,dti2,umol,dzz,grav,rho,kappa,u,v,dt,small,fsm,im,jm,kb,imm1,jmm1,kbm1,tbias,sbias,dz,...
         wusurf,wubot,wvsurf,wvbot,t,s,rhoref,zz,z)
% **********************************************************************
% *                                        Updated: Sep. 24, 2003      *
% * FUNCTION    :  Solves for q2 (twice the turbulent kinetic energy), *
% *                q2l (q2 x turbulent length scale), km (vertical     *
% *                kinematic viscosity) and kh (vertical kinematic     *
% *                diffusivity), using a simplified version of the     *
% *                level 2 1/2 model of Mellor and Yamada (1982).      *
% * In this version, the Craig-Banner sub-model whereby breaking wave  * 
% * tke is injected into the surface is included. However, we use an   *
% * analytical solution to the near surface tke equation to solve for  *
% * q2 at the surface giving the same result as C-B diffusion. The new *
% * scheme is simpler and more robust than the latter scheme.          *     
% *                                                                    *
% * References                                                         *
% *   Craig, P. D. and M. L. Banner, Modeling wave-enhanced turbulence *
% *     in the ocean surface layer. J. Phys. Oceanogr., 24, 2546-2559, *
% *     1994.                                                          *
% *   Ezer, T., On the seasonal mixed-layer simulated by a basin-scale *
% *     ocean model and the Mellor-Yamada turbulence scheme,           *
% *     J. Geophys. Res., 105(C7), 16,843-16,855, 2000.                *
% *   Mellor, G.L. and T. Yamada, Development of a turbulence          *
% *     closure model for geophysical fluid fluid problems,            *
% *     Rev. Geophys. Space Phys., 20, 851-875, 1982.                  *
% *   Mellor, G. L., One-dimensional, ocean surface layer modeling,    *
% *     a problem and a solution. J. Phys. Oceanogr., 31(3), 790-809,  *
% *     2001.                                                          *
% *   Mellor, G.L. and A. Blumberg, Wave breaking and ocean surface    *
% *     thermal response, J. Phys. Oceanogr., 2003.                    *
% *   Stacey, M. W., Simulations of the wind-forced near-surface       *
% *     circulation in Knight Inlet: a parameterization of the         *
% *     roughness length. J. Phys. Oceanogr., 29, 1363-1367, 1999.     *
% *                                                                    *
% **********************************************************************
%load('grid.mat');load('operator.mat');load('para.mat');
global dzz_3d kbm2 dz_3d zz_3d z_3d;

a1=0.92; b1=16.6 ; a2=0.74 ; b2=10.1    ; c1=0.08;
e1=1.8 ; e2=1.33 ; sef=1.0 ; cbcnst=100.; surfl=2.e5 ; shiw=0.0;
      
l0=zeros(im,jm);
kn = zeros(im,jm,kb);
boygr=zeros(im,jm,kb);
gh=zeros(im,jm,kb);
stf=zeros(im,jm,kb);
la=zeros(kb);
      
dh = h + etf;
dh_3d=repmat(dh,1,1,kb);

a = zeros(im,jm,kb);
c = zeros(im,jm,kb);
d = zeros(im,jm,kb);

a(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,2:kbm1);
c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,1:kbm2);

a= DIVISION(-dti2 .* (AZF(kq)+umol) , (a .* dh_3d .* dh_3d));
c= DIVISION(-dti2 .* (AZB(kq)+umol) , (c .* dh_3d .* dh_3d));
a(:,:,1)=0.e0;         c(:,:,1)=0.e0;       
a(:,:,kb)=0.e0;        c(:,:,kb)=0.e0;
%-----------------------------------------------------------------------
%     The following section solves the equation:
%
%       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
%
%     Surface and bottom boundary conditions:
%
const1=(16.6e0^(2.e0/3.e0))*sef;

% initialize fields that are not calculated on all boundaries
% but are later used there
l0(:,jm)  =0; l0(im,:)=0;
kn(:,:,:)=0;
      
utau2 = sqrt( AXF(wusurf).^2 +AYF(wvsurf).^2 );
% Wave breaking energy- a variant of Craig & Banner (1994), see Mellor and Blumberg, 2003.
 % Surface length scale following Stacey (1999).
l0 = surfl*utau2/grav;                     
uf(:,:,kb) = sqrt( AXF(wubot).^2 +AYF(wvbot).^2 ) .* const1;
        
%    Calculate speed of sound squared:
h_3d=repmat(h,1,1,kb);      
p=grav*rhoref*(-zz_3d .* h_3d)*1.e-4;     %     Calculate pressure in units of decibars:
cc=1449.10+.00821*p+4.55*(t+tbias) -.045e0*(t+tbias).^2 +1.34*(s+sbias-35.0e0);
cc=cc./sqrt((1.0-.01642.*p./cc) .*(1.0-0.40.*p./cc.^2));      
cc(:,:,kb)=0;
      
%     Calculate buoyancy gradient:
q2b =abs(q2b);
q2lb=abs(q2lb);
tmp = zeros(im,jm,kb);
for k=2:kbm1
    tmp(:,:,k)=dzz(k-1);
end
boygr=DIVISION(-grav* DZB(rho) , tmp .* h_3d ) + DIVISION(grav^2 , AZB(cc.^2));
boygr(:,:,1)=0.e0; boygr(:,:,kb)=0.e0;
l=q2lb ./ q2b;
%l=max(l, repmat(kappa*l0,1,1,kb));
l=(z_3d>-0.5) .* max(l, repmat(kappa*l0,1,1,kb))+(z_3d<=-0.5) .* l;
gh=l.^2 .* boygr ./q2b;
gh=min(gh, 0.028);

l(:,:,1)=kappa*l0; l(:,:,kb)=0;
gh(:,:,1)=0      ; gh(:,:,kb)=0;

%    Calculate production of turbulent kinetic energy:
kn= DIVISION(km.*sef.*(AXF(DZB(u)).^2 + AYF(DZB(v)).^2) , (tmp.*dh_3d).^2) -shiw.*km.*boygr + kh.*boygr;
%
%  NOTE: Richardson # dep. dissipation correction (Mellor: 2001; Ezer, 2000),
%  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
ghc=-6.0e0;
stf=ones(im,jm,kb);
% It is unclear yet if diss. corr. is needed when surf. waves are included.
%           if(gh(i,j,k).lt.0.e0)
%    ...        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
%           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
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
uf(filter) = small;
vf(filter) = 0.1 * dt_3d(filter) * small;

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
       
      fsm_3d = repmat(fsm, 1, 1, kb);
      km(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
      kh(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
      km(:,1,:)=km(:,2,:).*fsm_3d(:,1,:);
      kh(:,1,:)=kh(:,2,:).*fsm_3d(:,1,:);
      
      km(im,:,:)=km(imm1,:,:).*fsm_3d(im,:,:);
      kh(im,:,:)=kh(imm1,:,:).*fsm_3d(im,:,:);
      km(1,:,:)=km(2,:,:).*fsm_3d(1,:,:);
      kh(1,:,:)=kh(2,:,:).*fsm_3d(1,:,:);
      
 end
 
