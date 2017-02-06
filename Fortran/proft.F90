subroutine proft(f, wfsurf, dti2)
  use dm
  use dm_op
  use grid 
  use input

  implicit none 
  integer          :: ierr

  real(kind=8), intent(in) :: dti2
  type(Matrix), intent(inout) :: f
  type(Matrix), intent(in) ::  wfsurf

  type(Matrix) :: dh_3d, swrad_3d, rad, LAA, LBB, f1 , tmp_d 

  ! function [f1] = new_proft(f,wfsurf,fsurf,nbc,h,etf,swrad,kh)
  ! %-----------------------------------------------------------------------
  ! %     Irradiance parameters after Paulson and Simpson:
  ! %       ntp               1      2       3       4       5
  ! %   Jerlov type           i      ia      ib      ii     iii
  ! r=[0.58,0.62,0.67,0.77,0.78];
  ! ad1=[0.35,0.60,1.0,1.5,1.4];
  ! ad2=[0.23,20.0,17.0,14.0,7.90];
  ! % solves the equation:dti2*(kh*f')'-f=-fb

  ! rad=zeros(im,jm,kb);
  ! a = zeros(im,jm,kb);c = zeros(im,jm,kb); d=zeros(im,jm,kb);
  ! dh = h+etf;
  ! dh_3d=repmat(dh,1,1,kb);swrad_3d=repmat(swrad,1,1,kb);
  ! la=zeros(kbm1); f1=zeros(im,jm,kb);
  rad      = dm_zeros(im, jm, kb)
  dh_3d    = dm_rep(h + etf, 1, 1, kb)
  swrad_3d = dm_rep(swrad, 1, 1, kb)

  ! d(:,:,1:kbm2)=-dti2*(kh(:,:,2:kbm1)+umol);
  ! a=DIVISION(d,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
  ! a(:,:,kbm1)= 0.e0;
  tmp_d = -dti2 * ((shift(kh, 3, 1) + umol) .em. shift(REV_MASK_Z2, 3, 1))
  a = tmp_d .ed. (dz_3d .em. dzz_3d .em. dh_3d .em. dh_3d)

  ! d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
  ! c=DIVISION(-dti2*(kh+umol),dz_3d.*d.*dh_3d.*dh_3d);
  ! c(:,:,1)  =0.e0;
  tmp_d = tmp_d .em. MASK_Z1 + shift(dzz_3d, 3, -1) .em. REV_MASK_Z2
  c = (-dti2 * (kh + umol)) .ed. (dz_3d .em. tmp_d .em. dh_3d .em. dh_3d) &
       .em. REV_MASK_Z1

  !     d=-f - dti2 .* DZF(rad) ./(dh_3d .* dz_3d);
  !     d(:,:,1)= -f(:,:,1) + dti2 .* wfsurf(:,:) ./ (dh(:,:) .* dz(1));
  !     d(:,:,kb-1)=-f(:,:,kb-1); 
  tmp_d = -f - dti2 * DZF(rad) .ed. (dh_3d .em. dz_3d)
  tmp_d = tmp_d .em. REV_MASK_Z1 + &
       (-f + dti2 * dm_rep(wfsurf, 1, 1, kb) .ed. &
       (dh_3d .em.dz_3d)) .em. MASK_Z1

  tmp_d = tmp_d .em. (1 - shift(MASK_Z2, 3, 1)) + &
       (-f .em. shift(MASK_Z2, 3, 1))

  !   for j=2:jm
  !       for i=2:im
  !    la=diag(reshape(a(i,j,1:kbm1)+c(i,j,1:kbm1)-1,kbm1,1),0) ...
  !       - diag(reshape(a(i,j,1:kbm2),kbm2,1),1) ...
  !       - diag(reshape(c(i,j,2:kbm1),kbm2,1),-1);
  !    f1(i,j,1:kbm1)=la\reshape(d(i,j,1:kbm1),kbm1,1); 
  !       end
  !   end
  LAA = dm_trid(c, a + c - 1, a)
  LBB = dm_trid1(tmp_d)
  f1 = dm_trid2(dm_solve(LAA, LBB), im, jm, kb) 

  call dm_destroy(dh_3d, ierr)
  call dm_destroy(swrad_3d, ierr)
  call dm_destroy(rad, ierr)
  call dm_destroy(LAA, ierr)
  call dm_destroy(LBB, ierr)
  call dm_destroy(tmp_d, ierr)
end subroutine
