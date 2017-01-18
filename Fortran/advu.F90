subroutine advu(dti2)
  use dm
  use dm_op
  use grid
  use input
    implicit none
    integer        :: ierr
    type(Matrix)   :: tmpuf
    real(kind=8),intent(in) :: dti2

    h_3d     = dm_rep(h, 1, 1, kb)
    dt_3d    = dm_rep(dt, 1, 1, kb)
    etb_3d   = dm_rep(etb, 1, 1, kb)
    cor_3d   = dm_rep(cor, 1, 1, kb)
    egf_3d   = dm_rep(egf, 1, 1, kb)
    egb_3d   = dm_rep(egb, 1, 1, kb)
    e_atmos_3d = dm_rep(e_atmos, 1, 1, kb)
    etf_3d   = dm_rep(etf, 1, 1, kb)

    ! tmpuf       = dm_zeros(im, jm, kb)
    tmpuf = AXB(w) .em. AZB(u)

    uf = (AXB(etb_3d + h_3d) .em. aru_3d .em. ub - &
         dti2 * (advx+drhox-aru_3d .em. AXB(cor_3d .em.&
         dt_3d .em. AYF(v)) + 0.5e0*grav*AXB(dt_3d) .em.&
         (DXB(egf_3d+egb_3d)+2.e0*DXB(e_atmos_3d)) .em.&
         AXB(dy_3d) - DZF(tmpuf) .em. aru_3d .ed.&
         dz_3d)) .ed. (AXB(etf_3d + h_3d) .em. aru_3d) 
    uf = uf .em. REV_MASK_X2 + tmpuf .em. MASK_X2 
    call dm_destroy(tmpuf, ierr)    
  end subroutine 

