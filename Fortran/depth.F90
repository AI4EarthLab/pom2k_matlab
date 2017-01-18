  subroutine depth(ierr)
    use dm_op
    use dm
    use dm_type
    use input
    use grid
    implicit none

    integer, intent(out) :: ierr
    real(kind=8) :: kdz(12), delz
    integer :: i,k
    real(kind=8), allocatable :: z1(:), zz1(:), dz1(:), dzz1(:)

    allocate(z1(kb), zz1(kb), dz1(kb), dzz1(kb))

    kdz=(/1,1,2,4,8,16,32,64,128,256,512,1024/);

    z1 = 0; zz1 = 0; dz1 = 0; dzz1 = 0;

    do k=2,kl1
       z1(k)=z1(k-1)+kdz(k-1);
    enddo

    delz=z1(kl1)-z1(kl1-1);

    do k=kl1+1,kl2
       z1(k)=z1(k-1)+delz;
    enddo

    do k=kl2+1,kb
       dz1(k)=kdz(kb-k+1)*delz/kdz(kb-kl2);
       z1(k)=z1(k-1)+dz1(k);
    enddo

    do k=1,kb
       z1(k)=-z1(k)/z1(kb);
    enddo

    do k=1,kb-1
       zz1(k)=0.5e0*(z1(k)+z1(k+1));
    enddo

    zz1(kb)=2.e0*zz1(kb-1)-zz1(kb-2);

    do k=1,kb-1
       dz1(k)=z1(k)-z1(k+1);
       dzz1(k)=zz1(k)-zz1(k+1);
    enddo

    dz1(kb)=0.e0;
    dzz1(kb)=0.e0;

    ! print*, "z=", z1
    ! print*, "zz=", zz1
    ! print*, "dz=", dz1
    ! print*, "ddz=", dzz1  

    z  = dm_zeros(1, 1, kb)
    zz = dm_zeros(1, 1, kb)
    dz = dm_zeros(1, 1, kb)
    dzz= dm_zeros(1, 1, kb)


    call dm_setvalues(z, (/0/), (/0/), (/(i, i=0,kb-1)/), z1, ierr)
    call dm_setvalues(zz, (/0/), (/0/), (/(i, i=0,kb-1)/), zz1, ierr)
    call dm_setvalues(dz, (/0/), (/0/), (/(i, i=0,kb-1)/), dz1, ierr)
    call dm_setvalues(dzz, (/0/), (/0/), (/(i, i=0,kb-1)/), dzz1, ierr)


    z_3d  = dm_rep(z,  im, jm, 1)
    zz_3d = dm_rep(zz, im, jm, 1)
    dz_3d = dm_rep(dz, im, jm, 1)
    dzz_3d = dm_rep(dzz, im, jm, 1)

    deallocate(z1, zz1, dz1, dzz1)

  end subroutine 
