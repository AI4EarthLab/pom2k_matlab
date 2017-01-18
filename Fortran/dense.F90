  subroutine dense()
    use dm_op
    use dm
    use grid
    use input

    implicit none

    integer :: ierr
    type(Matrix) :: rhor, tr, sr, tr2, tr3, tr4, p, cr

    tr = t + tbias;
    sr = s + sbias;
    tr2 = tr .em. tr;
    tr3 = tr2 .em. tr;
    tr4 = tr3 .em. tr;

    !% Approximate pressure in units of bars:
    !p = (-1.e-5 * grav * rhoref) * (zz_3d .em. h_3d);

    p = zz_3d .em. h_3d


    rhor= -0.157406e0 + 6.793952e-2 * tr - 9.095290e-3 * tr2 + 1.001685e-4 .em. tr3 &
         -1.120083e-6 .em. tr4 + 6.536332e-9 .em. tr4 .em. tr;

    rhor=rhor + (0.824493e0-4.0899e-3*tr+7.6438e-5*tr2-8.2467e-7*tr3 &
         +5.3875e-9*tr4) .em. sr+(-5.72466e-3+1.0227e-4*tr &
         -1.6546e-6*tr2) .em. dm_pow(dm_abs(sr), 1.5)

    ! rhor=rhor + (0.824493e0-4.0899e-3*tr+7.6438e-5*tr2-8.2467e-7*tr3 &
    !      +5.3875e-9*tr4) .em. sr+(-5.72466e-3+1.0227e-4*tr &
    !      -1.6546e-6*tr2) .em. dm_pow(dm_abs(sr), 1.5) + 4.8314e-4 * sr .em. sr;
    cr=1449.1e0+.0821e0 * p+4.55e0 * tr-.045e0*tr2+1.34e0*(sr-35.e0);

    rhor=rhor +1.e5 * p .ed. (cr .em. cr) .em. (1.e0-2.e0 * p .ed. (cr .em. cr));

    rho= (rhor * (1.0/rhoref)) .em. fsm_3d;

    call dm_destroy(rhor, ierr)
    call dm_destroy(tr, ierr)
    call dm_destroy(sr, ierr)
    call dm_destroy(tr2, ierr)
    call dm_destroy(tr3, ierr)
    call dm_destroy(tr4, ierr)
    call dm_destroy(p, ierr)
    call dm_destroy(cr, ierr)

  end subroutine dense
