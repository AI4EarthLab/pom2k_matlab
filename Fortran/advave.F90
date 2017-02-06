  subroutine advave(iint)
    use dm_op
    use grid
    use input
    implicit none
    integer      :: ierr
    integer, intent(in) :: iint

    ! if(iint == 3) then
    !    call dm_print_info(d, ierr)
    !    call dm_print_info(aam2d, ierr)
    !    call dm_print_info(uab, ierr)
    !    call dm_print_info(vab, ierr)              
    !    call dm_print_info(dx, ierr)
    !    call dm_print_info(dy, ierr)
    !    call dm_finalize(ierr)
    !    stop
    !  endif

    print *, " in file: ", __FILE__, "line:", __LINE__

       call dm_print_info(d, ierr)
       call dm_print_info(aam2d, ierr)
       call dm_print_info(uab, ierr)
       call dm_print_info(vab, ierr)              
       call dm_print_info(dx, ierr)
       call dm_print_info(dy, ierr)
    
    tps   = AYB(AXB(d)) .em. AXB(AYB(aam2d)) .em. &
         (DYB(uab) .em. dm_pow(AYB(AXB(dy)), -1) + &
         DXB(vab) .em. dm_pow(AYB(AXB(dx)), -1))

    print *, " in file: ", __FILE__, "line:", __LINE__
    
    advua = DXB((AXF(AXB(d) .em. ua) .em. AXF(ua) - 2.e0*d .em. aam2d &
         .em. DXF(uab) .ed. dx) .em. dy) + &
         DYF((AXB(AYB(d) .em. va) .em. AYB(ua) - tps) .em. &
         AYB(AXB(dx)))

    advva = DYB((AYF(AYB(d) .em. va) .em. AYF(va) - 2.e0*d .em. aam2d &
         .em. DYF(vab) .ed. dy) .em. dx) + &
         DXF((AYB(AXB(d) .em. ua) .em. AXB(va) - tps) .em. &
         AYB(AXB(dy)))
    print *, " in file: ", __FILE__, "line:", __LINE__
    
    if (mode == 2) then
       print *, " in file: ", __FILE__, "line:", __LINE__
       
       wubot = -AXB(cbc) .em. dm_sqrt(dm_squ(uab)+dm_squ(AXB(AYF(vab)))) &
            .em. uab
       wvbot = -AYB(cbc) .em. dm_sqrt(dm_squ(vab)+dm_squ(AYB(AXF(uab)))) &
            .em. vab
       curv2d= (AYF(va) .em. DXC(dy) - AXF(ua) .em. DYC(dx)) .ed. &
            (dx .em. dy)
       advua = advua - aru .em. AXB(curv2d .em. d .em. AYF(va))
       advva = advva + arv .em. AYB(curv2d .em. d .em. AXF(ua))
       
       print *, " in file: ", __FILE__, "line:", __LINE__           
    end if

  end subroutine 
