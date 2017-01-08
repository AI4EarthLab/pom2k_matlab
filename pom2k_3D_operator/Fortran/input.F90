module input
  implicit none

  integer :: iend, im, jm, kb, kl1, kl2
  real(kind=8) :: dte, isplit, sw
  integer :: mode, nadv, nitera, nread
  real(kind=8) :: prtd1, prtd2, swtch
  integer :: iskp, jskp
  logical :: lramp
  real(kind=8) :: tprni, slmax
  integer :: ntp, nbct, nbcs, ispadv
  real(kind=8) :: alpha, tbias, sbias, tatm, satm
  integer :: iproblem
  
  namelist /setting/ iend, im, jm, kb, dte, isplit, kl1, kl2, iproblem, mode, nadv 
  namelist /setting/ nitera, sw, nread, prtd1, prtd2
  namelist /setting/ swtch, iskp, jskp, lramp, tprni, slmax, ntp 
  namelist /setting/ nbct, nbcs , ispadv, alpha , tbias, sbias, tatm, satm

contains
  
  subroutine LoadInput()
    open(1,file='input')
    read(1, setting)
    print*, "im=", im
  end subroutine  
  
end module input







