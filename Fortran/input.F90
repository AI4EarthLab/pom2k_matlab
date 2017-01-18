module input
  implicit none

  integer :: iend, im, jm, kb, kl1, kl2
  integer :: kbm1, kbm2, imm1, jmm1, imm2, jmm2
  real(kind=8) :: dte, isplit, sw
  integer :: mode, nadv, nitera, nread
  real(kind=8) :: prtd1, prtd2, swtch
  integer :: iskp, jskp
  logical :: lramp
  real(kind=8) :: tprni, slmax
  integer :: ntp, nbct, nbcs, ispadv
  real(kind=8) :: alpha, tbias, sbias, tatm, satm
  integer :: iproblem
  integer :: days
  real(kind=8) :: small, pi, rhoref, grav, kappa
  real(kind=8) :: z0b, cbcmin, cbcmax, horcon, umol, hmax
  real(kind=8) :: smoth, aam_init, ramp, vmaxl
  
  namelist /setting/ iend, im, jm, kb, dte, isplit, kl1, kl2, iproblem, mode, nadv 
  namelist /setting/ nitera, sw, nread, prtd1, prtd2, days
  namelist /setting/ swtch, iskp, jskp, lramp, tprni, slmax, ntp
  namelist /setting/ nbct, nbcs , ispadv, alpha , tbias, sbias, tatm, satm
  namelist /setting/ small, pi, rhoref, grav, kappa, z0b
  namelist /setting/ cbcmin, cbcmax, horcon, umol, hmax, vmaxl
  namelist /setting/ smoth, aam_init, ramp
  
contains
  
  subroutine LoadInput()
    open(1, file='input')
    read(1, setting)
    imm1 = im - 1
    imm2 = im - 2
    jmm1 = jm - 1
    jmm2 = jm - 2
    kbm1 = kb - 1
    kbm2 = kb - 2
  end subroutine  
  
end module input







