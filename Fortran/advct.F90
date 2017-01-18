subroutine advct()
  use dm
  use dm_op
  use grid
  use input
  implicit none

  integer         :: ierr

  ! curv = dm_zeros(im, jm, kb)
  ! xflux= dm_zeros(im, jm, kb)
  ! yflux= dm_zeros(im, jm, kb)

  ! advx = dm_zeros(im, jm, kb)
  ! advy = dm_zeros(im, jm, kb)

  curv = (AYF(v) .em. DXB(AXF(dy_3d)) - &
       AXF(u) .em. DYB(AYF(dx_3d))) .ed. (dx_3d .em. dy_3d)

  xflux= dy_3d .em. (AXF(AXB(dt_3d) .em. u) .em. AXF(u) &
       - 2.e0 * dt_3d .em. aam .em. DXF(ub) .ed. dx_3d)

  yflux= AYB(AXB(dx_3d)) .em. ((AXB(AYB(dt_3d) .em. v) .em. AYB(u)) &
       - AYB(AXB(dt_3d)) .em. AYB(AXB(aam)) .em. &
       (DYB(ub) .ed. AYB(AXB(dy_3d)) + DXB(vb) .ed. AYB(AXB(dx_3d))))

  advx = DXB(xflux) + DYF(yflux) - aru_3d .em. AXB(curv .em. dt_3d .em. AYF(v))

  yflux= dx_3d .em. (AYF(AYB(dt_3d) .em. v) .em. AYF(v) &
       - 2.e0 * dt_3d .em. aam .em. DYF(vb) .ed. dy_3d)

  xflux= AYB(AXB(dy_3d)) .em. ((AYB(AXB(dt_3d) .em. u) .em. AXB(v)) &
       - AYB(AXB(dt_3d)) .em. AYB(AXB(aam)) .em. &
       (DYB(ub) .ed. AYB(AXB(dy_3d)) + DXB(vb) .ed. AYB(AXB(dx_3d))))

  advy = DXF(xflux) + DYB(yflux) + arv_3d .em. &
       AYB(curv .em. dt_3d .em. AXF(u)) ! Modify "-aru_3d" to "+arv_3d""

end subroutine advct
