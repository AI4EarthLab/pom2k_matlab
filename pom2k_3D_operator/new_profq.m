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
load('grid.mat');load('operator.mat');load('para.mat');

a1=0.92; b1=16.6 ; a2=0.74 ; b2=10.1    ; c1=0.08;
e1=1.8 ; e2=1.33 ; sef=1.0 ; cbcnst=100.; surfl=2.e5 ; shiw=0.0;
      
l0=zeros(im,jm);
kn = zeros(im,jm,kb);
boygr=zeros(im,jm,kb);
gh=zeros(im,jm,kb);
stf=zeros(im,jm,kb);
      
dh = h + etf;
dh_3d=repmat(dh,1,1,kb);

a = zeros(im,jm,kb);
c = zeros(im,jm,kb);
tmp1 = zeros(im,jm,kb);
tmp2 = zeros(im,jm,kb);
for k=2:kbm1
    tmp1(:,:,k)=dzz(k-1)*dz(k);
    tmp2(:,:,k)=dzz(k-1)*dz(k-1);
end

a= DIVISION(-dti2 .* (AZF1(kq)+umol) , (tmp1 .* dh_3d .* dh_3d));
c= DIVISION(-dti2 .* (AZB1(kq)+umol) , (tmp2 .* dh_3d .* dh_3d));
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
ee(:,jm,1)=0; ee(im,:,1)=0; 
gg(:,jm,1)=0; gg(im,:,1)=0; 
l0(:,jm)  =0; l0(im,:)=0;
kn(:,:,:)=0;
      
utau2 = sqrt( AXF1(wusurf).^2 +AYF1(wvsurf).^2 );
% Wave breaking energy- a variant of Craig & Banner (1994), see Mellor and Blumberg, 2003.
ee(:,:,1)=zeros(im,jm);
gg(:,:,1)=(15.8*cbcnst)^(2./3.)*utau2;
 % Surface length scale following Stacey (1999).
l0 = surfl*utau2/grav;                     
uf(:,:,kb) = sqrt( AXF1(wubot).^2 +AYF1(wvbot).^2 ) .* const1;
        
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
boygr=DIVISION(-grav* DZB1(rho) , tmp .* h_3d ) + DIVISION(grav^2 , AZB2(cc.^2));

l=q2lb ./ q2b;
l=max(l, repmat(kappa*l0,1,1,kb));
gh=l.^2 .* boygr ./q2b;
gh=min(gh, 0.028);

l(:,:,1)=kappa*l0; l(:,:,kb)=0;
gh(:,:,1)=0      ; gh(:,:,kb)=0;

%    Calculate production of turbulent kinetic energy:
kn= DIVISION(km.*sef.*(AXF2(DZB2(u)).^2 + AYF2(DZB2(v)).^2) , (tmp.*dh_3d).^2) -shiw.*km.*boygr + kh.*boygr
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

%
for k=2:kbm1
    gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-(2.e0*dti2.*dtef(:,:,k)+1.e0));
    ee(:,:,k)=a(:,:,k).*gg(:,:,k);
    gg(:,:,k)=(-2.e0*dti2*kn(:,:,k)+c(:,:,k).*gg(:,:,k-1)-uf(:,:,k)).*gg(:,:,k);
end

for k=1:kbm1
    ki=kb-k;
    uf(:,:,ki)=ee(:,:,ki).*uf(:,:,ki+1)+gg(:,:,ki);
end      
%-----------------------------------------------------------------------
%     The following section solves the equation:
%
%       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
ee(:,:,2)=0.e0;
gg(:,:,2)=0.e0;
vf(:,:,kb)=0.e0;

for k=2:kbm1
    dtef(:,:,k)=dtef(:,:,k).*(1.e0+e2.*((1.e0/abs(z(k)-z(1))+1.e0./abs(z(k)-z(kb)))     ...
                        .*l(:,:,k)./(dh(:,:)*kappa)).^2);
    gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-(dti2*dtef(:,:,k)+1.e0));
    ee(:,:,k)=a(:,:,k).*gg(:,:,k);
    gg(:,:,k)=(dti2*(-kn(:,:,k).*l(:,:,k).*e1)+c(:,:,k).*gg(:,:,k-1)-vf(:,:,k)).*gg(:,:,k);
end

for k=1:kb-2
    ki=kb-k;
    vf(:,:,ki)=ee(:,:,ki).*vf(:,:,ki+1)+gg(:,:,ki);
end

for k=2:kbm1
    for j=1:jm
        for i=1:im
            if(uf(i,j,k)<=small || vf(i,j,k)<=small) 
                uf(i,j,k)=small;
                vf(i,j,k)=0.1*dt(i,j)*small;
            end
        end
    end
end
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
% cosmetics: make boundr. values as interior
% (even if not used: printout otherwise may show strange values)
      for k=1:kb
        for i=1:im
           km(i,jm,k)=km(i,jmm1,k)*fsm(i,jm);
           kh(i,jm,k)=kh(i,jmm1,k)*fsm(i,jm);
           km(i,1,k)=km(i,2,k)*fsm(i,1);
           kh(i,1,k)=kh(i,2,k)*fsm(i,1);
        end
        for j=1:jm
           km(im,j,k)=km(imm1,j,k)*fsm(im,j);
           kh(im,j,k)=kh(imm1,j,k)*fsm(im,j);
           km(1,j,k)=km(2,j,k)*fsm(1,j);
           kh(1,j,k)=kh(2,j,k)*fsm(1,j);
        end
      end
 end