module constants

implicit none
integer, parameter :: im        = 65
integer, parameter :: jm        = 49
integer, parameter :: kb        = 21
real(8), parameter :: small     = 1.d-9                                       ! a small value
real(8), parameter :: pi        = 3.141592653589793238462643383279502884      ! PI
real(8), parameter :: rhoref    = 1025.0                                      ! Reference density (recommended values: 1025 for
                                                                              ! seawater, 1000 for freshwater, S.I. units)
real(8), parameter :: grav      = 9.8060                                      ! gravity constant (S.I. units)
real(8), parameter :: kappa     = 0.40                                        ! von Karman's constant'
real(8), parameter :: z0b       = 0.01                                        ! Bottom roughness (metres)
real(8), parameter :: cbcmin    = 0.0025                                      ! Minimum bottom friction coeff
real(8), parameter :: cbcmax    = 1.0                                         ! Maximum bottom friction coeff
real(8), parameter :: horcon    = 0.2                                         ! Smagorinsky diffusivity coeff
real(8), parameter :: umol      = 2.d-5                                       ! Background viscosity used in subroutines profq,
                                                                              ! proft, profu, profv (S.I. units)
real(8), parameter :: hmax      = 4500.0                                      ! Maximum depth used in radiation boundary condition 
                                                                              ! in subroutine bcond(meters)
real(8), parameter :: vmax      = 100.0                                       ! Maximum magnitude of vaf (used in check that
                                                                              ! essentially tests for CFL violation)
real(8), parameter :: smoth     = 0.10                                        ! Constant in temporal filter used to prevent
                                                                              ! solution splitting (dimensionless)
real(8), parameter :: aam_init  = 500.0                                       ! Initial value of aam.
real(8), parameter :: ramp      = 0.0                                         ! A coeff which is not used in original Fortran.

end module constants







