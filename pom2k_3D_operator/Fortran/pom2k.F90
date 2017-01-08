program pom2k
  use dm_op 
  use dm
  use dm_type
  use constants
  use input
  use grid
  implicit none
  type(Matrix) :: tmp
  integer :: ierr
  integer :: myrank, mysize

  call dm_init1(ierr)
  call dm_comm_rank(myrank,ierr)
  call dm_comm_size(mysize,ierr)

  call InitGrid()
  call InitOperatorModule(im, jm, kb)

  call dm_init2(im, jm , kb,ierr)

  if(myrank==0) then
     print *, "==============Input paramenters==========="
     print *, "im=",im,",jm=",jm,",kb=",kb
  endif


  call new_depth(z, zz, dz, dzz, z_3d, zz_3d, dz_3d, dzz_3d, kl1, kl2, ierr)

  call new_dense()
  
  call new_baropg()

  !print*, tmp%nx, tmp%ny, tmp%nz
  !call dm_view( z_3d .em. MASK_Z1, ierr)

  ! call dm_view(z, ierr)
  ! call dm_view(zz, ierr)
  ! call dm_view(dz, ierr)
  ! call dm_view(dzz, ierr)

  ! call dm_view(z_3d, ierr)
  ! call dm_view(zz_3d, ierr)
  ! call dm_view(dz_3d, ierr)
  ! call dm_view(dzz_3d, ierr)

  call FinalizeOperatorModule()

  call FinalizeGrid()  
  call dm_finalize(ierr)

end program pom2k
