subroutine advv(dti2)
  use dm
  use dm_op
  use grid
  use input
  implicit none
  integer        :: ierr
  real(kind=8),intent(in) :: dti2
  type(Matrix) :: tmpvf

  h_3d     = dm_rep(h, 1, 1, kb)
  dt_3d    = dm_rep(dt, 1, 1, kb)
  etb_3d   = dm_rep(etb, 1, 1, kb)
  cor_3d   = dm_rep(cor, 1, 1, kb)
  egf_3d   = dm_rep(egf, 1, 1, kb)
  egb_3d   = dm_rep(egb, 1, 1, kb)
  e_atmos_3d = dm_rep(e_atmos, 1, 1, kb)
  etf_3d   = dm_rep(etf, 1, 1, kb)

  ! vf       = dm_zeros(im, jm, kb)
  tmpvf = AYB(w) .em. AZB(v)

  vf = (AYB(etb_3d + h_3d) .em. arv_3d .em. vb - &
       dti2 * (advy+drhoy+arv_3d .em. AYB(cor_3d .em.&
       dt_3d .em. AXF(u)) + 0.5e0*grav*AYB(dt_3d) .em.&
       (DYB(egf_3d+egb_3d)+2.e0*DYB(e_atmos_3d)) .em.&
       AYB(dx_3d) - DZF(tmpvf) .em. arv_3d .ed.&
       dz_3d)) .ed. (AYB(etf_3d + h_3d) .em. arv_3d) 

  call dm_destroy(tmpvf, ierr)

end subroutine advv
