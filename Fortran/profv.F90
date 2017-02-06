subroutine profv(dti2)
  use dm
  use dm_op
  use grid
  use input

  implicit none
  real(kind=8) :: dti2
  type(Matrix) ::  LAA, LBB, dh, dh_3d, cbc_3d, tps_3d, wvsurf_3d, tmp_d
  integer :: i, ierr
  !  a=zeros(im,jm,kb);
  ! dh = AYB(h+etf); dh(1,:)=1.e0; dh(:,1)=1.e0;
  dh = AYB(h + etf)
  call dm_setvalues(dh, (/0/), (/(i,i=0,jm-1)/), (/0/), (/(1, i=0, jm-1)/), ierr)
  call dm_setvalues(dh, (/(i,i=0,im-1)/), (/0/), (/0/), (/(1, i=0, im-1)/), ierr)

  print *,  " in file: ", __FILE__, "line:", __LINE__
  
  ! dh_3d=repmat(dh,1,1,kb);
  dh_3d = dm_rep(dh, 1, 1, kb)
  
  ! c = AYB(km); c(1,:,:)=0.e0;c(:,1,:)=0.e0;
  c = AYB(km) .em. REV_MASK_X1 .em. REV_MASK_Y1
  
  ! la=zeros(kbm1);d=zeros(im,jm,kb);
  ! %
  !     d(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
  !     a=DIVISION(d,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
  !     a(:,:,kbm1)= 0.e0;
  tmp_d = -dti2 * shift(c+umol, 3, 1)
  a = tmp_d .ed. (dz_3d .em. dzz_3d .em. dh_3d .em. dh_3d)
  a = a .em. (1-shift(MASK_Z2, 3, 1))
  
  !     d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
  !     c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);
  !     c(:,:,1)  =0.e0;
  tmp_d = (shift(dzz_3d, 3, 1) .em. MASK_Z2) + (tmp_d .em. (MASK_ZZ))
  c = (-dti2 * (c+umol)) .ed. (dz_3d .em. tmp_d .em. dh_3d .em. dh_3d)
  c = c .em. REV_MASK_Z1
  
  ! tps = AYB(cbc) .* sqrt( AYB( AXF( ub(:,:,kbm1) ) ).^2 + vb(:,:,kbm1).^2 );
  cbc_3d = dm_rep(cbc, 1, 1, kb)
  tps_3d = AYB(cbc_3d) .em. dm_sqrt(dm_squ(AYB(AXF(ub))) + dm_squ(vb)) .em. MASK_Z2
  
  !     d=-vf;
  !     d(:,:,1)= -vf(:,:,1) + dti2 .* wvsurf(:,:) ./ (dh(:,:) .* dz(1));
  !     d(:,:,kbm1)=-vf(:,:,kbm1)+tps.*dti2./(dz(kbm1)*dh); 
  tmp_d = -vf
  wvsurf_3d = dm_rep(wvsurf, 1, 1, kb)
  dh_3d = dm_rep(dh, 1, 1, kb)
  
  tmp_d = tmp_d .em. REV_MASK_Z1 + &
       (MASK_Z1 .em. (dti2 * (wvsurf_3d .ed. (dh_3d .em. dz_3d))))
  
  tmp_d = tmp_d .em. (1-shift(MASK_Z1,3,-1)) + &
       shift(MASK_Z2, 3, -1) .em. (dti2 * (tps_3d .ed. (dh_3d .em. dz_3d)))
  
  !   for j=2:jm
  !       for i=2:im
  !    la=diag(reshape(a(i,j,1:kbm1)+c(i,j,1:kbm1)-1,kbm1,1),0) ...
  !       - diag(reshape(a(i,j,1:kbm2),kbm2,1),1) ...
  !       - diag(reshape(c(i,j,2:kbm1),kbm2,1),-1);
  !    vf(i,j,1:kbm1)=la\reshape(d(i,j,1:kbm1),kbm1,1); 
  !       end
  !   end
  !    vf=vf.*repmat(dvm,1,1,kb);
  !    wvbot=-tps .* vf(:,:,kbm1);
  LAA = dm_trid(c, a + c - 1, a)
  LBB = dm_trid1(tmp_d)
  vf  = dm_trid2(dm_solve(LAA, LBB), im, jm, kb) .em. dm_rep(dvm, 1, 1, kb)
  wvbot = csum(shift(MASK_Z2, 3, -1) .em. (-tps_3d) .em. vf, 4)
  
  call dm_destroy(LAA, ierr)
  call dm_destroy(LBB, ierr)
  call dm_destroy(dh, ierr)
  call dm_destroy(dh_3d, ierr)
  call dm_destroy(cbc_3d, ierr)
  call dm_destroy(tps_3d, ierr)
  call dm_destroy(wvsurf_3d, ierr)
  call dm_destroy(tmp_d, ierr)
end subroutine profv
