module constants

implicit none
! integer, parameter :: im        = 65
! integer, parameter :: jm        = 49
! integer, parameter :: kb        = 21
real(8), parameter :: ! a small value
real(8), parameter :: 653589793238462643383279502884      ! PI
real(8), parameter :: ! Reference density (recommended values: 1025 for
                      ! seawater, 1000 for freshwater, S.I. units)
real(8), parameter :: ! gravity constant (S.I. units)
real(8), parameter :: ! von Karman's constant'
real(8), parameter :: ! Bottom roughness (metres)
real(8), parameter :: ! Minimum bottom friction coeff
real(8), parameter :: ! Maximum bottom friction coeff
real(8), parameter :: ! Smagorinsky diffusivity coeff
real(8), parameter ::  Background viscosity used in subroutines profq,
                      ! proft, profu, profv (S.I. units)
real(8), parameter :: ! Maximum depth used in radiation boundary condition 
                      ! in subroutine bcond(meters)
real(8), parameter :: ! Maximum magnitude of vaf (used in check that
                      ! essentially tests for CFL violation)
real(8), parameter :: ! Constant in temporal filter used to prevent
                      ! solution splitting (dimensionless)
real(8), parameter :: ! Initial value of aam.
real(8), parameter :: ! A coeff which is not used in original Fortran.

end module constants







