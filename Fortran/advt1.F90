!function advt1(fb, f, dti2) result(ff)
subroutine advt1(ff, fb, f, dti2)
  use dm
  use dm_op
  use grid
  use input
  implicit none
  integer  :: ierr
  type(Matrix), intent(in)  :: fb, f
  type(Matrix), intent(out) :: ff
  real(kind=8), intent(in)  :: dti2

  dt_3d = dm_rep(dt, 1, 1, kb)
  h_3d  = dm_rep(h, 1, 1, kb)

  !xflux=(AXB(dt_3d).*AXB(f).*u-DIVISION(AXB(aam).*  ...
  !       AXB(h_3d)*tprni.*DXB(fb).*dum_3d , AXB(dx_3d))).*AXB(dy_3d);
  ! yflux=(AYB(dt_3d).*AYB(f).*v-DIVISION(AYB(aam).* ...
  !        AYB(h_3d)*tprni.*DYB(fb).*dvm_3d , AYB(dy_3d))).*AYB(dx_3d);
  xflux = (AXB(dt_3d) .em. AXB(f) .em. u - (tprni * AXB(aam) .em.&
       AXB(h_3d) .em. DXB(fb) .em. dum_3d) .ed. AXB(dx_3d)) .em.&
       AXB(dy_3d)
  yflux = (AYB(dt_3d) .em. AYB(f) .em. v - (tprni * AYB(aam) .em.&
       AYB(h_3d) .em. DYB(fb) .em. dvm_3d) .ed. AYB(dy_3d)) .em.&
       AYB(dx_3d)

  zflux = AZB(f) .em. w .em. art_3d

  !xflux(:, :, kb) = 0.e0
  !yflux(:, :, kb) = 0.e0
  !zflux(:, :, kb) = 0.e0
  xflux = xflux .em. REV_MASK_Z2
  yflux = yflux .em. REV_MASK_Z2
  zflux = zflux .em. REV_MASK_Z2

  ff = (fb .em. (h_3d + dm_rep(etb, 1, 1, kb)) .em. art_3d - &
       dti2 * (DXF(xflux) + DYF(yflux) - DZF(zflux) .ed. dz_3d))&
       .ed. ((h_3d + dm_rep(etf, 1, 1, kb)) .em. art_3d)

  !ff(1, :, :) = 0.e0
  !ff(im, :, :)= 0.e0
  !ff(:, 1, :) = 0.e0
  !ff(:, jm, :)= 0.e0
  !ff(:, :, kb)= 0.e0
  ff = ff .em. REV_MASK_X1
  ff = ff .em. REV_MASK_X2
  ff = ff .em. REV_MASK_Y1
  ff = ff .em. REV_MASK_Y2
  ff = ff .em. REV_MASK_Z2
end subroutine
