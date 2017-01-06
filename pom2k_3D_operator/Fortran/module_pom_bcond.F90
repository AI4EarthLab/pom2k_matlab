module pom_bcond

use dm

implicit none

contains

subroutine bcond1(elf, ierr)        ! For external (2-D) boundary condition, elevation.
    implicit none
    type(Matrix), intent(inout) :: elf

end subroutine bcond1
end module pom_cond
