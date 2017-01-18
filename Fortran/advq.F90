subroutine advq(qf, qb, q, dti2) 
  use dm
  use dm_op
  use grid
  use input
  implicit none
  type(Matrix), intent(inout) :: qb, q, qf
  real(kind=8), intent(in) :: dti2

  integer        :: ierr

  !type(Matrix)   :: xflux, yflux
  ! xflux = dm_zeros(im, jm, kb)
  ! yflux = dm_zeros(im, jm, kb)
  ! qf    = dm_zeros(im, jm, kb)

  dt_3d = dm_rep(dt, 1, 1, kb)
  h_3d  = dm_rep(h, 1, 1, kb)

  xflux = AXB(dy_3d) .em. (AXB(q) .em. AXB(dt_3d) .em. AZB(u) -&
       AZB(AXB(aam)) .em. AXB(h_3d) .em. DXB(qb) .em.&
       dum_3d .em. dm_pow(AXB(dx_3d), -1))

  !call dm_print_info(AXB(aam), ierr)
  !call dm_view(AXB(aam), ierr)
  !xflux = AZB(AXB(aam))
  ! xflux = AXB(dy_3d) .em. (AXB(q) .em. AXB(dt_3d) .em. AZB(u) - &
  !       AXB(AXB(aam))) !.em. AXB(h_3d) .em. DXB(qb) .em.&
  ! dum_3d .em. dm_pow(AXB(dx_3d), -1))

  yflux = AYB(dx_3d) .em. (AYB(q) .em. AYB(dt_3d) .em. AZB(v) -&
       AZB(AYB(aam)) .em. AYB(h_3d) .em. DYB(qb) .em.&
       dvm_3d .em. dm_pow(AYB(dy_3d), -1))

  qf    = ((h_3d + dm_rep(etb, 1, 1, kb)) .em. art_3d .em. qb - &
       dti2 * (-DZC(w .em. q) .em. art_3d .em. &
       dm_pow(AZB(dz_3d), -1) + DXF(xflux)+DYF(yflux))) .ed. &
       ((h_3d + dm_rep(etf, 1, 1, kb)) .em. art_3d)

  ! qf(1, :, :) = 0.e0
  ! qf(im, :, :)= 0.e0
  ! qf(:, 1, :) = 0.e0
  ! qf(:, jm, :)= 0.e0
  ! qf(:, :, 1) = 0.e0
  ! qf(:, :, kb)= 0.e0      
  qf = qf .em. REV_MASK_X1
  qf = qf .em. REV_MASK_X2
  qf = qf .em. REV_MASK_Y1
  qf = qf .em. REV_MASK_Y2
  qf = qf .em. REV_MASK_Z1
  qf = qf .em. REV_MASK_Z2

end subroutine

