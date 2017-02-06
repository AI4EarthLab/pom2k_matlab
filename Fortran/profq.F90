subroutine profq(dti2)
  use dm
  use dm_op
  use grid
  use input
  implicit none

  ! function [sm,sh,dh,cc,ee,gg,l,kq,km,kh,...
  !     uf,vf,q2b,q2lb,a,c]=sx_profq(sm,sh,dh,cc,...
  !     ee,gg,l,kq,km,kh,uf,vf,q2,q2b,q2lb,a,c,...
  !     h,etf,dti2,umol,dzz,grav,rho,kappa,u,v,dt,small,fsm,
  !     im,jm,kb,imm1,jmm1,kbm1,tbias,sbias,dz,...
  !     wusurf,wubot,wvsurf,wvbot,t,s,rhoref,zz,z)

  real(kind=8), intent(in) :: dti2
  real(kind=8) ::  a1=0.92, b1=16.6, a2=0.74, b2=10.1, c1=0.08, const1
  real(kind=8) ::  e1=1.8, e2=1.33, sef=1.0, cbcnst=100., surfl=2.e5, shiw=0.0
  type(Matrix) ::  l0, kn, boygr, gh, dh_3d2, dh, dh_3d,cc,p,utau2
  type(Matrix) :: filter, l0_kappa, tmp
  real(kind=8) :: ghc, coef1, coef2, coef3, coef4, coef5
  type(Matrix) :: LAA, LBB, sh, sm, temp, temp1, temp2, tmp_d
  integer :: ierr

  !l0    = dm_zeros(im,jm,1);
  !kn    = dm_zeros(im,jm,kb);
  !boygr = dm_zeros(im,jm,kb);
  !gh    = dm_zeros(im,jm,kb);
  !stf   = dm_zeros(im,jm,kb);


  dh = h + etf;
  dh_3d=dm_rep(dh,1,1,kb);

  a = dm_zeros(im,jm,kb);
  c = dm_zeros(im,jm,kb);
  !d = dm_zeros(im,jm,kb);


  ! a(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,2:kbm1);
  ! c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,1:kbm2);
  a = SHIFT(dzz_3d, 3, -1) .em. dz_3d
  c = dzz_3d .em. dz_3d


  ! a= DIVISION(-dti2 .* (AZF(kq)+umol) , (a .* dh_3d .* dh_3d));
  ! c= DIVISION(-dti2 .* (AZB(kq)+umol) , (c .* dh_3d .* dh_3d));
  dh_3d2 = dm_squ(dh_3d)
  a = (-dti2 * (AZF(kq) + umol)) .ed. (a .em. dh_3d2)
  c = (-dti2 * (AZB(kq) + umol)) .ed. (c .em. dh_3d2)

  ! a(:,:,1)=0.e0;         c(:,:,1)=0.e0;       
  ! a(:,:,kb)=0.e0;        c(:,:,kb)=0.e0;
  a = a .em. REV_MASK_ZZ
  c = c .em. REV_MASK_ZZ

  ! %-----------------------------------------------------------------------
  ! %     The following section solves the equation:
  ! %       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
  ! %     Surface and bottom boundary conditions:
  const1=(16.6e0**(2.e0/3.e0))*sef;


  ! % initialize fields that are not calculated on all boundaries
  ! % but are later used there
  ! l0(:,jm)  =0; l0(im,:)=0;
  ! kn(:,:,:)=0;

  ! utau2 = sqrt( AXF(wusurf).^2 +AYF(wvsurf).^2 );
  utau2 = dm_sqrt(dm_squ(AXF(wusurf)) + dm_squ(AYF(wvsurf)));


  !  % Surface length scale following Stacey (1999).
  ! l0 = surfl*utau2/grav;
  l0 = (surfl/grav) * utau2



  ! uf(:,:,kb) = sqrt( AXF(wubot).^2 +AYF(wvbot).^2 ) .* const1;
  uf = uf .em. REV_MASK_Z2 &
       + dm_rep(const1 * dm_sqrt(dm_squ(AXF(wubot)) + dm_squ(AYF(wvbot))), 1, 1, kb) &
       .em. REV_MASK_Z1


  ! %    Calculate speed of sound squared:
  ! h_3d=repmat(h,1,1,kb);
  h_3d = dm_rep(h, 1,1,kb)

  !%     Calculate pressure in units of decibars:
  ! p=grav*rhoref*(-zz_3d .* h_3d)*1.e-4;
  p = -grav*rhoref*1.e-4*(zz_3d .em. h_3d)

  ! cc=1449.10+.00821*p+4.55*(t+tbias) -.045e0*(t+tbias).^2 +1.34*(s+sbias-35.0e0);
  cc = 1449.10+0.00821*p+4.55*(t+tbias)-.045e0*dm_squ(t+tbias) + &
       1.34*(s+sbias-35.0e0)

  ! cc=cc./sqrt((1.0-.01642.*p./cc) .*(1.0-0.40.*p./cc.^2));
  cc = cc .ed. dm_sqrt((1.0-0.01642*p .ed. cc) .em. (1.0-0.40*p .ed. dm_squ(cc)))
  ! cc(:,:,kb)=0;
  cc = cc .em. REV_MASK_Z2


  ! %     Calculate buoyancy gradient:
  ! q2b =abs(q2b);
  q2b = dm_abs(q2b)

  ! q2lb=abs(q2lb);
  q2lb = dm_abs(q2lb)

  ! tmp = zeros(im,jm,kb);
  ! for k=2:kbm1
  !     tmp(:,:,k)=dzz(k-1);
  ! end
  tmp = SHIFT(dzz_3d, 3, -1)

  ! boygr=DIVISION(-grav* DZB(rho) , tmp .* h_3d ) + DIVISION(grav^2 , AZB(cc.^2));
  boygr = (-grav*DZB(rho)) .ed. (tmp .em. h_3d) + grav**2 * dm_pow(AZB(dm_squ(cc)), -1)

  ! boygr(:,:,1)=0.e0; boygr(:,:,kb)=0.e0;
  boygr = boygr .em. REV_MASK_ZZ

  ! l=q2lb ./ q2b;
  l = q2lb .ed. q2b

  ! %l=max(l, repmat(kappa*l0,1,1,kb));
  ! l=(z_3d>-0.5) .* max(l, repmat(kappa*l0,1,1,kb))+(z_3d<=-0.5) .* l;
  l0_kappa = dm_rep(kappa * l0, 1, 1, kb)
  filter = (l > l0_kappa)
  l = (z_3d > -0.5) .em. (filter .em. l + (1 - filter) .em. l0_kappa) + &
       (z_3d <= -0.5) .em. l

  ! print *, "in file: ", __FILE__, "line:", __LINE__
  ! call dm_finalize(ierr)
  ! stop

  ! gh=l.^2 .* boygr ./q2b;
  ! gh=min(gh, 0.028);
  gh = dm_squ(l) .em. boygr .ed. q2b
  filter = (gh >= 0.028)
  gh = 0.028 * filter + gh .em. (1-filter)

  ! l(:,:,1)=kappa*l0; l(:,:,kb)=0;
  ! gh(:,:,1)=0      ; gh(:,:,kb)=0;
  l = l .em. REV_MASK_ZZ
  gh = gh .em. REV_MASK_ZZ

  ! %    Calculate production of turbulent kinetic energy:
  ! kn= DIVISION(km.*sef.*(AXF(DZB(u)).^2 + AYF(DZB(v)).^2) , ...
  !             (tmp.*dh_3d).^2) -shiw.*km.*boygr + kh.*boygr;
  kn = (sef * km .em. (dm_squ(AXF(DZB(u))) + dm_squ(AYF(DZB(v))))) .ed. &
       (dm_squ(tmp .em. dh_3d) - shiw * km .em. boygr + kh .em. boygr)


  ! %  NOTE: Richardson # dep. dissipation correction (Mellor: 2001; Ezer, 2000),
  ! %  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
  ghc=-6.0e0;
  !stf=dm_ones(im,jm,kb);

  print *, "in file: ", __FILE__, "line:", __LINE__
  
  ! % It is unclear yet if diss. corr. is needed when surf. waves are included.
  ! %           if(gh(i,j,k).lt.0.e0)
  ! %    ...        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
  ! %           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
  ! dtef=sqrt(q2b).*stf./(b1.*l+small);
  ! dtef(:,:,1)=0.e0;dtef(:,:,kb)=0.e0;
  dtef = dm_sqrt(q2b) .ed. (b1 * l + small) .em. REV_MASK_ZZ

  !     d=-uf-2.e0*dti2*kn;
  !     d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
  !     d(:,jm,1)=0; d(im,:,1)=0; d(:,:,kb)=-uf(:,:,kb);
  !     temp1=a(:,:,1:kbm1);    temp2=c(:,:,2:kb);
  tmp_d = -2.0*dti2*kn-uf
  tmp_d = tmp_d .em. REV_MASK_Z1 - (15.8*cbcnst)**(2./3.) * dm_rep(utau2, 1, 1, kb)
  tmp_d = (1-(MASK_Y2 .em. MASK_Z1)) .em. tmp_d
  tmp_d = (1-(MASK_X2 .em. MASK_Z1)) .em .tmp_d
  tmp_d = (REV_MASK_Z1 .em. tmp_d) - (uf .em. MASK_Z1)
  temp1 = a .em. REV_MASK_Z2
  temp2 = c .em. REV_MASK_Z1

  print *, "in file: ", __FILE__, "line:", __LINE__

  ! call dm_finalize(ierr)
  ! stop

  !   for j=2:jm
  !       for i=2:im
  !    la=diag(reshape(a(i,j,:)+c(i,j,:)-1-2.e0*dti2.*dtef(i,j,:),kb,1),0) ...
  !       - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
  !       - diag(reshape(temp2(i,j,:),kbm1,1),-1);
  !    uf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
  !       end
  !   end

  !(m,n,k) => (k, k, m*n) 
  LAA = dm_trid(temp2, a+c-1-2.0*dti2*dtef, temp1)
  !(m,n,k) => (k, 1, m*n)
  LBB = dm_trid1(tmp_d)
  uf = dm_trid2(dm_solve(LAA, LBB), im, jm, kb)

  ! %-----------------------------------------------------------------------
  ! %     The following section solves the equation:
  ! %
  ! %       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
  ! vf(:,:,kb)=0.e0;
  vf = vf .em. MASK_Z2

  print *, "in file: ", __FILE__, "line:", __LINE__
  
  ! for k=2:kbm1
  !     dtef(:,:,k)=dtef(:,:,k).*(1.e0+e2.*((1.e0/abs(z(k)-z(1))+1.e0./abs(z(k)-z(kb)))     ...
  !                         .*l(:,:,k)./(dh(:,:)*kappa)).^2);
  ! end
  !     d=-vf-dti2*kn.*l.*e1;
  !     d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
  !     d(:,jm,1)=0.e0; d(im,:,1)=0.e0; d(:,:,kb)=0.e0;
  !     temp=vf(:,:,1);

  ! call dm_print_info(kn, ierr)
  ! call dm_print_info(l, ierr)
  ! call dm_print_info(vf, ierr)
  ! call dm_print_info(d, ierr)
  ! call dm_print_info(REV_MASK_X1, ierr)
  
  tmp_d = -dti2 * e1 * kn .em. l - vf
  tmp_d = dm_rep(-(15.8*cbcnst)**(2./3.) * utau2, 1, 1, kb) .em. MASK_Z1 +&
       (REV_MASK_Z1 .em. tmp_d)
  tmp_d = tmp_d .em. REV_MASK_X1
  tmp_d = tmp_d .em. REV_MASK_Y1
  tmp_d = tmp_d .em. REV_MASK_Z1
  temp = vf .em. MASK_Z1
  
  print *, "in file: ", __FILE__, "line:", __LINE__
  
  !   for j=2:jm
  !       for i=2:im
  !    la=diag(reshape(a(i,j,:)+c(i,j,:)-1-dti2.*dtef(i,j,:),kb,1),0) ...
  !       - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
  !       - diag(reshape(temp2(i,j,:),kbm1,1),-1);
  !    vf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
  !       end
  !   end
  LAA = dm_trid(temp2, a+c-1-dti2 * dtef, temp1)
  print *, "in file: ", __FILE__, "line:", __LINE__  
  LBB = dm_trid1(tmp_d)
  print *, "in file: ", __FILE__, "line:", __LINE__  
  vf = dm_trid2(dm_solve(LAA, LBB), im, jm, kb)
  print *, "in file: ", __FILE__, "line:", __LINE__
  !    vf(:,:,1)=temp;
  vf = vf .em. REV_MASK_Z1 + temp

  ! dt_3d=repmat(dt,1,1,kb);
  dt_3d = dm_rep(dt, 1, 1, kb)

  print *, "in file: ", __FILE__, "line:", __LINE__
  
  ! filter = (uf<=small | vf<=small);
  ! filter(:,:,1) = false;
  ! filter(:,:,kb) = false;
  ! uf(filter) = small;
  ! vf(filter) = 0.1 * dt_3d(filter) * small;
  filter = 1-((1 - (uf <= small)) .em. (1 - (vf <= small)))
  filter = filter .em. REV_MASK_Z1
  filter = filter .em. REV_MASK_Z2
  uf = (uf .em. (1-filter)) + (small * filter)
  vf = 0.1 * small * (dt_3d .em. filter) + (vf .em. (1-filter))

  ! %-----------------------------------------------------------------------
  ! %     The following section solves for km and kh:
  !     coef4=18.e0*a1*a1+9.e0*a1*a2;
  !     coef5=9.e0*a1*a2;
  ! %     Note that sm and sh limit to infinity when gh approaches 0.0288:
  !     coef1=a2*(1.e0-6.e0*a1/b1*stf);
  !     coef2=3.e0*a2*b2/stf+18.e0*a1*a2;
  !     coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf);
  !     sh=coef1./(1.e0-coef2.*gh);
  !     sm=coef3+sh.*coef4.*gh;
  !     sm=sm./(1.0-coef5.*gh);
  coef4=18.e0*a1*a1+9.e0*a1*a2;
  coef5=9.e0*a1*a2;
  coef1=a2*(1.e0-6.e0*a1/b1);
  coef2=3.e0*a2*b2+18.e0*a1*a2;
  coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1);

  sh=coef1 * dm_pow(1.e0-coef2*gh, -1);
  sm=coef3 + coef4 * sh .em. gh;
  sm=sm .em. dm_pow(1.0-coef5 * gh, -1);

  print *, "in file: ", __FILE__, "line:", __LINE__
  
  !     kn=l.*sqrt(abs(q2));
  !     kq=(kn.*.41e0.*sh+kq)*.5e0;
  !     km=(kn.*sm+km)*.5e0;
  !     kh=(kn.*sh+kh)*.5e0;
  kn = l .em. dm_sqrt(dm_abs(q2))
  kq = (0.41*kn .em. sh + kq) * 0.5
  km = (kn .em. sm + km) * 0.5
  kh = (kn .em. sh + kh) * 0.5

  !       fsm_3d = repmat(fsm, 1, 1, kb);
  !       km(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
  !       kh(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
  !       km(:,1,:)=km(:,2,:).*fsm_3d(:,1,:);
  !       kh(:,1,:)=kh(:,2,:).*fsm_3d(:,1,:);
  !fsm_3d = dm_rep(fsm, 1, 1, kb)

  !       km(im,:,:)=km(imm1,:,:).*fsm_3d(im,:,:);
  !       kh(im,:,:)=kh(imm1,:,:).*fsm_3d(im,:,:);
  !       km(1,:,:)=km(2,:,:).*fsm_3d(1,:,:);
  !       kh(1,:,:)=kh(2,:,:).*fsm_3d(1,:,:);
  km = (km .em. REV_MASK_X2) + (MASK_X2 .em. SHIFT(km, 1, -1) .em. fsm_3d)        
  kh = (kh .em. REV_MASK_X2) + (MASK_X2 .em. SHIFT(kh, 1, -1) .em. fsm_3d)    
  km = km .em. REV_MASK_X1 + MASK_X1 .em. SHIFT(km, 1, 1) .em. fsm_3d    
  kh = kh .em. REV_MASK_X1 + MASK_X1 .em. SHIFT(kh, 1, 1) .em. fsm_3d
  !  end
  print *, "in file: ", __FILE__, "line:", __LINE__
  
  call dm_destroy(l0, ierr)
  call dm_destroy(kn, ierr)
  call dm_destroy(boygr, ierr)
  call dm_destroy(gh, ierr)

  call dm_destroy(dh_3d2, ierr)
  call dm_destroy(dh, ierr)
  call dm_destroy(dh_3d, ierr)
  call dm_destroy(cc, ierr)
  call dm_destroy(p, ierr)
  call dm_destroy(utau2, ierr)
  call dm_destroy(tmp_d, ierr)
  
  print *, "in file: ", __FILE__, "line:", __LINE__
end subroutine profq

