  function [kq,km,kh,uf,vf,q2b,q2lb]=new_profq(ee,gg,l,kq,km,kh,uf,vf,q2,q2b,q2lb,...
                                        h,etf,dzz,rho,u,v,dt,fsm,wusurf,wubot,wvsurf,wvbot,t,s)
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
%
      
%
%      real sm(im,jm,kb),sh(im,jm,kb),cc(im,jm,kb)
%      real gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
%      real prod(im,jm,kb),kn(im,jm,kb)
%      real a1,a2,b1,b2,c1
%      real coef1,coef2,coef3,coef4,coef5
%      real const1,e1,e2,ghc
%      real p,sef,sp,tp
%      real l0(im,jm)
%      real cbcnst,surfl,shiw
%      real utau2, df0,df1,df2 
 %     equivalence (prod,kn)
%
load('grid.mat');load('operator.mat');load('para.mat');

sm=zeros(im,jm,kb);
sh=zeros(im,jm,kb);
dh=zeros(im,jm,kb);

      a1=0.92e0;b1=16.6e0;a2=0.74e0;b2=10.1e0;c1=0.08e0;
      e1=1.8e0;e2=1.33e0;
      sef=1.e0;
      cbcnst=100.;surfl=2.e5;shiw=0.0;
      
% 
      l0=zeros(im,jm);
      kn = zeros(im,jm);
      boygr=zeros(im,jm,kb);
      gh=zeros(im,jm,kb);
      stf=zeros(im,jm,kb);
      
      dh = h + etf;
      
%%%%%%---------------------------------------------------      
a = zeros(im,jm,kb);
 c = zeros(im,jm,kb);
 temp1 = zeros(im,jm,kb);
 temp2 = zeros(im,jm,kb);
 for k=2:kbm1
     temp1(:,:,k)=dzz(k-1)*dz(k);
     temp2(:,:,k)=dzz(k-1)*dz(k-1);
 end
 for j=1:jm
     a(:,j,:) = -dti2*(AZF2_XZ(permute(kq(:,j,:),[1,3,2]))+umol)...
                ./(permute(temp1(:,j,:),[1,3,2]).*(repmat(dh(:,j),1,kb).^2));
     c(:,j,:) = -dti2*(AZB2_XZ(permute(kq(:,j,:),[1,3,2]))+umol)...
                ./(permute(temp2(:,j,:),[1,3,2]).*(repmat(dh(:,j),1,kb).^2));
 end
a(isinf(a))=0;
c(isinf(c))=0;

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
  
    
kn = zeros(im,jm,kb);
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

cc=zeros(im,jm,kb);
 for k=1:kbm1
%     Calculate pressure in units of decibars:
%
            p=grav*rhoref*(-zz(k)* h)*1.e-4;
            cc(:,:,k)=1449.1e0+.00821e0*p+4.55e0*tp(:,:,k) -.045e0*tp(:,:,k).^2     ...
                 +1.34e0*(sp(:,:,k)-35.0e0);
            cc(:,:,k)=cc(:,:,k)     ...
                 ./sqrt((1.e0-.01642e0*p./cc(:,:,k))     ...  
                 .*(1.e0-0.40e0*p./cc(:,:,k).^2));
 end
%
%     Calculate buoyancy gradient:
%
%
q2b=abs(q2b);
q2lb=abs(q2lb);
temp = zeros(im,jm,kb);
 for k=2:kbm1
     temp(:,:,k)=dzz(k-1);
 end
for j=1:jm
    boygr(:,j,:)=-grav*DZB1_XZ(permute(rho(:,j,:),[1,3,2]))   ...
                    ./(permute(temp(:,j,:),[1,3,2]).*repmat(h(:,j),1,kb));
% *** NOTE: comment out next line if dens fores not include pressure     ... 
     +(grav^2)./AZB1_XZ(permute(cc(:,j,:),[1,3,2]).^2); 
end   

%
%l=abs(q2lb./q2b);
      for k=2:kbm1
        for j=1:jm
          for i=1:im
             l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k));
            if(z(k)>-0.5)
                 l(i,j,k)=max(l(i,j,k),kappa*l0(i,j));
            end
            gh(i,j,k)=(l(i,j,k)^2)*boygr(i,j,k)/q2b(i,j,k);
            gh(i,j,k)=min(gh(i,j,k),.028e0);
          end
        end
      end
%

          l(:,:,1)=kappa*l0;
          l(:,:,kb)=0.e0;
          gh(:,:,1)=0.e0;
          gh(:,:,kb)=0.e0;

%
%    Calculate production of turbulent kinetic energy:
%

for k=2:kbm1
            kn(:,:,k)=km(:,:,k)*sef     ...  
                 .*(AXF2_XY( u(:,:,k)-u(:,:,k-1) ).^2     ...
                 +AYF2_XY( v(:,:,k)-v(:,:,k-1) ).^2)./( (dzz(k-1)*dh).^2 )
%   Add shear due to internal wave field     ...
             -shiw*km(:,:,k).*boygr(:,:,k);
            kn(:,:,k)=kn(:,:,k)+kh(:,:,k).*boygr(:,:,k);
            kn(:,:,k)=OP_L_XY*kn(:,:,k);
end
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
    
for k=2:kbm1
          gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))     ...
                      -(2.e0*dti2*dtef(:,:,k)+1.e0));
          ee(:,:,k)=a(:,:,k).*gg(:,:,k);
          gg(:,:,k)=(-2.e0*dti2*kn(:,:,k)+c(:,:,k).*gg(:,:,k-1)     ...
                 -uf(:,:,k)).*gg(:,:,k);
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
          ee(:,:,2)=0.e0;
          gg(:,:,2)=0.e0;
          vf(:,:,kb)=0.e0;
%
for k=2:kbm1
            dtef(:,:,k)=dtef(:,:,k)     ...
                   .*(1.e0+e2*((1.e0/abs(z(k)-z(1))     ...
                               +1.e0/abs(z(k)-z(kb)))     ...
                                .*l(:,:,k)./(dh*kappa)).^2);
            gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))     ...
                      -(dti2*dtef(:,:,k)+1.e0));
            ee(:,:,k)=a(:,:,k) .* gg(:,:,k);
            gg(:,:,k)=(dti2*(-kn(:,:,k) .* l(:,:,k)*e1)     ...
                 +c(:,:,k).*gg(:,:,k-1)-vf(:,:,k)).*gg(:,:,k);
end
%
for k=1:kb-2
        ki=kb-k;
            vf(:,:,ki)=ee(:,:,ki).*vf(:,:,ki+1)+gg(:,:,ki);
end
%
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
 end
