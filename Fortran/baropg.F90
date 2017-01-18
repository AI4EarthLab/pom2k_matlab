subroutine baropg()
  use dm
  use dm_op
  use grid
  use input
  implicit none
  integer :: ierr

  ! new_baropg(rho, rmean, dt, ramp)
  ! % **********************************************************************
  ! % *                                                                    *
  ! % * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
  ! % *                                                                    *
  ! % **********************************************************************
  ! global im jm kb zz_3d dx_3d dy_3d grav dum_3d dvm_3d;

  rho = rho - rmean;
  drhox=dm_zeros(im,jm,kb);
  drhoy=dm_zeros(im,jm,kb);

  dt_3d=dm_rep(dt,1,1,kb);

  ! !% Calculate x-component of baroclinic pressure gradient:
  ! ! drhox(:,:,1)= -grav .* zz_3d(:,:,1) .* AXB(dt_3d(:,:,1)) .* DXB(rho(:,:,1));  
  drhox = drhox + (NAG_MASK_Z1 .em. drhox) - &
       grav .em. (zz_3d .em. MASK_Z1) .em. &
       AXB(dt_3d .em. MASK_Z1) .em. DXB(rho .em. MASK_Z1)

  ! drhox= SUM1(drhox-grav *DZB(zz_3d) .* AXB(dt_3d) .* DXB(AZB(rho))
  !        + grav * AZB(zz_3d) .* DXB(dt_3d) .* DZB(AXB(rho)));
  drhox = CSUM(drhox - grav .em. DZB(zz_3d) .em. AXB(dt_3d) .em. DXB(AZB(rho)) &
       + grav .em. AZB(zz_3d) .em. DXB(dt_3d) .em. DZB(AXB(rho)), 1)

  ! drhox = ramp*drhox .* AXB(dt_3d) .* dum_3d .* AXB(dy_3d);
  drhox = ramp * drhox .em. AXB(dt_3d) .em. dum_3d .em. AXB(dy_3d);


  ! %Calculate y-component of baroclinic pressure gradient:
  ! drhoy(:,:,1)= -grav .* zz_3d(:,:,1) .* AYB(dt_3d(:,:,1)) .* DYB(rho(:,:,1));
  drhoy = drhoy + (drhoy .em. NAG_MASK_Z1) - &
       grav .em. (zz_3d .em. MASK_Z1) .em. &
       AYB(dt_3d .em. MASK_Z1) .em. DYB(rho .em. MASK_Z1)

  ! drhoy= SUM1(drhoy-grav * DZB(zz_3d) .* AYB(dt_3d) .*
  ! DYB(AZB(rho)) + grav * AZB(zz_3d) .* DYB(dt_3d) .* DZB(AYB(rho)));
  drhoy = CSUM(drhoy - grav .em. DZB(zz_3d) .em. AYB(dt_3d) .em. DYB(AZB(rho)) &
       + grav * AZB(zz_3d) .em. DYB(dt_3d) .em. DZB(AYB(rho)), 1)

  ! drhoy = ramp*drhoy .* AYB(dt_3d) .* dvm_3d .* AYB(dx_3d);
  drhoy = ramp * drhoy .em. AYB(dt_3d) .em. dvm_3d .em. AYB(dx_3d)

  ! rho = rho + rmean;
  rho = rho + rmean

end subroutine baropg
