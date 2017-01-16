module functions
  use dm_op
  use grid
  use input
  use dm

contains

  subroutine new_dense()
    use dm_op
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

  end subroutine new_dense

  function new_advq(qb, q, dti2) result(qf)
    implicit none
    type(Matrix), intent(inout) :: qb, q
    real(kind=8), intent(in) :: dti2
    type(Matrix) :: qf
    integer        :: ierr

    !type(Matrix)   :: xflux, yflux

    ! xflux = dm_zeros(im, jm, kb)
    ! yflux = dm_zeros(im, jm, kb)

    ! qf    = dm_zeros(im, jm, kb)

    ! call dm_print_info(dy_3d, ierr)
    ! call dm_print_info(q, ierr)
    ! call dm_print_info(dt_3d, ierr)
    ! call dm_print_info(u, ierr)
    ! call dm_print_info(aam, ierr)
    ! call dm_print_info(h_3d, ierr)
    ! call dm_print_info(qb, ierr)
    ! call dm_print_info(dum_3d, ierr)
    ! call dm_print_info(dx_3d, ierr)

    dt_3d = dm_rep(dt, 1, 1, kb)
    h_3d  = dm_rep(h, 1, 1, kb)

    xflux = AXB(dy_3d) .em. (AXB(q) .em. AXB(dt_3d) .em. AZB(u) -&
         AZB(AXB(aam)) .em. AXB(h_3d) .em. DXB(qb) .em.&
         dum_3d .em. dm_pow(AXB(dx_3d), -1))

    !call dm_print_info(AXB(aam), ierr)
    !call dm_view(AXB(aam), ierr)
    !xflux = AZB(AXB(aam))
    ! xflux = AXB(dy_3d) .em. (AXB(q) .em. AXB(dt_3d) .em. AZB(u) - &
    !       AXB(AXB(aam))) !.em. AXB(h_3d) .em. DXB(qb) .em.&
    ! dum_3d .em. dm_pow(AXB(dx_3d), -1))


    yflux = AYB(dx_3d) .em. (AYB(q) .em. AYB(dt_3d) .em. AZB(v) -&
         AZB(AYB(aam)) .em. AYB(h_3d) .em. DYB(qb) .em.&
         dvm_3d .em. dm_pow(AYB(dy_3d), -1))

    qf    = ((h_3d + dm_rep(etb, 1, 1, kb)) .em. art_3d .em. qb - &
         dti2 * (-DZC(w .em. q) .em. art_3d .em. &
         dm_pow(AZB(dz_3d), -1) + DXF(xflux)+DYF(yflux))) .ed. &
         ((h_3d + dm_rep(etf, 1, 1, kb)) .em. art_3d)


    ! qf(1, :, :) = 0.e0
    ! qf(im, :, :)= 0.e0
    ! qf(:, 1, :) = 0.e0
    ! qf(:, jm, :)= 0.e0
    ! qf(:, :, 1) = 0.e0
    ! qf(:, :, kb)= 0.e0      
    qf = qf .em. REV_MASK_X1
    qf = qf .em. REV_MASK_X2
    qf = qf .em. REV_MASK_Y1
    qf = qf .em. REV_MASK_Y2
    qf = qf .em. REV_MASK_Z1
    qf = qf .em. REV_MASK_Z2

    ! call dm_destroy(xflux, ierr)
    ! call dm_destroy(yflux, ierr)

  end function new_advq

  subroutine new_baropg()
    implicit none
    integer :: ierr

    ! new_baropg(rho, rmean, dt, ramp)
    ! % **********************************************************************
    ! % *                                                                    *
    ! % * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
    ! % *                                                                    *
    ! % **********************************************************************
    ! global im jm kb zz_3d dx_3d dy_3d grav dum_3d dvm_3d;

    rho = rho - rmean;
    drhox=dm_zeros(im,jm,kb);
    drhoy=dm_zeros(im,jm,kb);

    dt_3d=dm_rep(dt,1,1,kb);

    ! !% Calculate x-component of baroclinic pressure gradient:
    ! ! drhox(:,:,1)= -grav .* zz_3d(:,:,1) .* AXB(dt_3d(:,:,1)) .* DXB(rho(:,:,1));  
    drhox = drhox + (NAG_MASK_Z1 .em. drhox) - &
         grav .em. (zz_3d .em. MASK_Z1) .em. &
         AXB(dt_3d .em. MASK_Z1) .em. DXB(rho .em. MASK_Z1)

    ! drhox= SUM1(drhox-grav *DZB(zz_3d) .* AXB(dt_3d) .* DXB(AZB(rho))
    !        + grav * AZB(zz_3d) .* DXB(dt_3d) .* DZB(AXB(rho)));
    drhox = CSUM(drhox - grav .em. DZB(zz_3d) .em. AXB(dt_3d) .em. DXB(AZB(rho)) &
         + grav .em. AZB(zz_3d) .em. DXB(dt_3d) .em. DZB(AXB(rho)), 1)

    ! drhox = ramp*drhox .* AXB(dt_3d) .* dum_3d .* AXB(dy_3d);
    drhox = ramp * drhox .em. AXB(dt_3d) .em. dum_3d .em. AXB(dy_3d);


    ! %Calculate y-component of baroclinic pressure gradient:
    ! drhoy(:,:,1)= -grav .* zz_3d(:,:,1) .* AYB(dt_3d(:,:,1)) .* DYB(rho(:,:,1));
    drhoy = drhoy + (drhoy .em. NAG_MASK_Z1) - &
         grav .em. (zz_3d .em. MASK_Z1) .em. &
         AYB(dt_3d .em. MASK_Z1) .em. DYB(rho .em. MASK_Z1)

    ! drhoy= SUM1(drhoy-grav * DZB(zz_3d) .* AYB(dt_3d) .*
    ! DYB(AZB(rho)) + grav * AZB(zz_3d) .* DYB(dt_3d) .* DZB(AYB(rho)));
    drhoy = CSUM(drhoy - grav .em. DZB(zz_3d) .em. AYB(dt_3d) .em. DYB(AZB(rho)) &
         + grav * AZB(zz_3d) .em. DYB(dt_3d) .em. DZB(AYB(rho)), 1)

    ! drhoy = ramp*drhoy .* AYB(dt_3d) .* dvm_3d .* AYB(dx_3d);
    drhoy = ramp * drhoy .em. AYB(dt_3d) .em. dvm_3d .em. AYB(dx_3d)

    ! rho = rho + rmean;
    rho = rho + rmean

  end subroutine new_baropg

  subroutine new_advave(iint)

    implicit none

    integer      :: ierr
    integer, intent(in) :: iint
    !type(Matrix) :: tps

    ! advua = dm_zeros(im, jm, 1)
    ! advva = dm_zeros(im, jm, 1)
    ! tps   = dm_zeros(im, jm, 1)

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

    tps   = AYB(AXB(d)) .em. AXB(AYB(aam2d)) .em. &
         (DYB(uab) .em. dm_pow(AYB(AXB(dy)), -1) + &
         DXB(vab) .em. dm_pow(AYB(AXB(dx)), -1))

    
    advua = DXB((AXF(AXB(d) .em. ua) .em. AXF(ua) - 2.e0*d .em. aam2d &
         .em. DXF(uab) .ed. dx) .em. dy) + &
         DYF((AXB(AYB(d) .em. va) .em. AYB(ua) - tps) .em. &
         AYB(AXB(dx)))

    advva = DYB((AYF(AYB(d) .em. va) .em. AYF(va) - 2.e0*d .em. aam2d &
         .em. DYF(vab) .ed. dy) .em. dx) + &
         DXF((AYB(AXB(d) .em. ua) .em. AXB(va) - tps) .em. &
         AYB(AXB(dy)))



     
    if (mode == 2) then

       wubot = -AXB(cbc) .em. dm_sqrt(dm_squ(uab)+dm_squ(AXB(AYF(vab)))) &
            .em. uab
       wvbot = -AYB(cbc) .em. dm_sqrt(dm_squ(vab)+dm_squ(AYB(AXF(uab)))) &
            .em. vab
       curv2d= (AYF(va) .em. DXC(dy) - AXF(ua) .em. DYC(dx)) .ed. &
            (dx .em. dy)
       advua = advua - aru .em. AXB(curv2d .em. d .em. AYF(va))
       advva = advva + arv .em. AYB(curv2d .em. d .em. AXF(ua))

    end if

    !call dm_destroy(tps, ierr)
  end subroutine new_advave

  subroutine new_advct()
    implicit none

    integer         :: ierr

    ! curv = dm_zeros(im, jm, kb)
    ! xflux= dm_zeros(im, jm, kb)
    ! yflux= dm_zeros(im, jm, kb)

    ! advx = dm_zeros(im, jm, kb)
    ! advy = dm_zeros(im, jm, kb)

    curv = (AYF(v) .em. DXB(AXF(dy_3d)) - &
         AXF(u) .em. DYB(AYF(dx_3d))) .ed. (dx_3d .em. dy_3d)

    xflux= dy_3d .em. (AXF(AXB(dt_3d) .em. u) .em. AXF(u) &
         - 2.e0 * dt_3d .em. aam .em. DXF(ub) .ed. dx_3d)

    yflux= AYB(AXB(dx_3d)) .em. ((AXB(AYB(dt_3d) .em. v) .em. AYB(u)) &
         - AYB(AXB(dt_3d)) .em. AYB(AXB(aam)) .em. &
         (DYB(ub) .ed. AYB(AXB(dy_3d)) + DXB(vb) .ed. AYB(AXB(dx_3d))))

    advx = DXB(xflux) + DYF(yflux) - aru_3d .em. AXB(curv .em. dt_3d .em. AYF(v))

    yflux= dx_3d .em. (AYF(AYB(dt_3d) .em. v) .em. AYF(v) &
         - 2.e0 * dt_3d .em. aam .em. DYF(vb) .ed. dy_3d)

    xflux= AYB(AXB(dy_3d)) .em. ((AYB(AXB(dt_3d) .em. u) .em. AXB(v)) &
         - AYB(AXB(dt_3d)) .em. AYB(AXB(aam)) .em. &
         (DYB(ub) .ed. AYB(AXB(dy_3d)) + DXB(vb) .ed. AYB(AXB(dx_3d))))

    advy = DXF(xflux) + DYB(yflux) + arv_3d .em. &
         AYB(curv .em. dt_3d .em. AXF(u)) ! Modify "-aru_3d" to "+arv_3d""

  end subroutine new_advct

  subroutine new_vertvl(dti2)
    implicit none
    real(kind=8),intent(in) :: dti2
    integer :: ierr

    ! function [xflux,yflux,w]=
    ! new_vertvl (xflux,yflux,w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,
    !             etb,dti2,im,jm,imm1,jmm1,kbm1)
    ! % **********************************************************************
    ! % *                                                                    *
    ! % * FUNCTION    :  calculates vertical velocity.                       *
    ! % *                                                                    *
    ! % **********************************************************************

    ! dt_3d=repmat(dt,1,1,kb);
    dt_3d=dm_rep(dt,1,1,kb);

    !% Reestablish boundary conditions:
    ! xflux = AXB(dy_3d) .* AXB(dt_3d).* u;
    ! yflux = AYB(dx_3d) .* AYB(dt_3d).* v;
    a = AXB(dy_3d) .em. AXB(dt_3d) .em. u;
    c = AYB(dx_3d) .em. AYB(dt_3d) .em. v;

    ! %     NOTE that, if one wishes to include freshwater flux, the
    ! %     surface velocity should be set to vflux(i,j). See also
    ! %     change made to 2-D volume conservation equation which
    ! %     calculates elf.
    ! w(:,:,1)=0.5e0*(vfluxb+vfluxf);
    ! w=repmat(w(:,:,1),1,1,kb);              
    w  = dm_rep(0.5 * (vfluxb + vfluxf), 1, 1, kb)

    ! tps=(etf-etb)/dti2;
    ! tps=repmat(tps,1,1,kb);
    tps = (etf-etb) * (1.0/dti2);
    tps = dm_rep(tps,1,1,kb);


    !%by hx
    ! w=w+SUM2(dz_3d .* ( ( DXF(xflux)*1.e0+DYF(yflux) )./(dx_3d .*dy_3d)+ tps));  %by hx
    w = w+CSUM(dz_3d .em. ((DXF(a)*1.e0+DYF(c)) .ed. (dx_3d .em. dy_3d)+ tps), 2);  

    ! w(1,:,:) = 0.e0;
    ! w(im,:,:) = 0.e0;
    w = w .em. MASK_X1
    w = w .em. MASK_X2

  end subroutine new_vertvl

  function new_advt1(fb, f, dti2) result(ff)
    implicit none
    integer         :: ierr
    type(Matrix), intent(in) :: fb, f
    type(Matrix) :: ff
    real(kind=8), intent(in) :: dti2

    ! xflux = dm_zeros(im, jm, kb)
    ! yflux = dm_zeros(im, jm, kb)
    ! zflux = dm_zeros(im, jm, kb)

    ! ff    = dm_zeros(im, jm, kb)
    dt_3d = dm_rep(dt, 1, 1, kb)
    h_3d  = dm_rep(h, 1, 1, kb)

    xflux = (AXB(dt_3d) .em. AXB(f) .em. u - (tprni * AXB(aam) .em.&
         AXB(h_3d) .em. DXB(fb) .em. dum_3d) .ed. AXB(dx_3d)) .em.&
         AXB(dy_3d)
    yflux = (AYB(dt_3d) .em. AYB(f) .em. v - (tprni * AYB(aam) .em.&
         AYB(h_3d) .em. DYB(fb) .em. dvm_3d) .ed. AYB(dy_3d)) .em.&
         AYB(dx_3d)

    zflux = AZB(f) .em. w .em. art_3d

    !xflux(:, :, kb) = 0.e0
    !yflux(:, :, kb) = 0.e0
    !zflux(:, :, kb) = 0.e0
    xflux = xflux .em. REV_MASK_Z2
    yflux = yflux .em. REV_MASK_Z2
    zflux = zflux .em. REV_MASK_Z2

    ff = (fb .em. (h_3d + dm_rep(etb, 1, 1, kb)) .em. art_3d - &
         dti2 * (DXF(xflux) + DYF(yflux) - DZF(zflux) .ed. dz_3d))&
         .ed. ((h_3d + dm_rep(etf, 1, 1, kb)) .em. art_3d)

    !ff(1, :, :) = 0.e0
    !ff(im, :, :)= 0.e0
    !ff(:, 1, :) = 0.e0
    !ff(:, jm, :)= 0.e0
    !ff(:, :, kb)= 0.e0
    ff = ff .em. REV_MASK_X1
    ff = ff .em. REV_MASK_X2
    ff = ff .em. REV_MASK_Y1
    ff = ff .em. REV_MASK_Y2
    ff = ff .em. REV_MASK_Z2

  end function new_advt1

  subroutine new_advu(dti2)
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
  end subroutine new_advu

  subroutine new_advv(dti2)
    implicit none
    integer        :: ierr
    real(kind=8),intent(in) :: dti2
    type(Matrix) :: tmpvf

    h_3d     = dm_rep(h, 1, 1, kb)
    dt_3d    = dm_rep(dt, 1, 1, kb)
    etb_3d   = dm_rep(etb, 1, 1, kb)
    cor_3d   = dm_rep(cor, 1, 1, kb)
    egf_3d   = dm_rep(egf, 1, 1, kb)
    egb_3d   = dm_rep(egb, 1, 1, kb)
    e_atmos_3d = dm_rep(e_atmos, 1, 1, kb)
    etf_3d   = dm_rep(etf, 1, 1, kb)

    ! vf       = dm_zeros(im, jm, kb)
    tmpvf = AYB(w) .em. AZB(v)

    vf = (AYB(etb_3d + h_3d) .em. arv_3d .em. vb - &
         dti2 * (advy+drhoy+arv_3d .em. AYB(cor_3d .em.&
         dt_3d .em. AXF(u)) + 0.5e0*grav*AYB(dt_3d) .em.&
         (DYB(egf_3d+egb_3d)+2.e0*DYB(e_atmos_3d)) .em.&
         AYB(dx_3d) - DZF(tmpvf) .em. arv_3d .ed.&
         dz_3d)) .ed. (AYB(etf_3d + h_3d) .em. arv_3d) 

    call dm_destroy(tmpvf, ierr)

  end subroutine new_advv


  subroutine new_depth(ierr)
    use dm_op
    use dm
    use dm_type
    use input
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

  end subroutine new_depth

  subroutine new_profq(dti2)
    use dm
    use dm_op
    
    implicit none
    
    ! function [sm,sh,dh,cc,ee,gg,l,kq,km,kh,...
    !     uf,vf,q2b,q2lb,a,c]=sx_profq(sm,sh,dh,cc,...
    !     ee,gg,l,kq,km,kh,uf,vf,q2,q2b,q2lb,a,c,...
    !     h,etf,dti2,umol,dzz,grav,rho,kappa,u,v,dt,small,fsm,
    !     im,jm,kb,imm1,jmm1,kbm1,tbias,sbias,dz,...
    !     wusurf,wubot,wvsurf,wvbot,t,s,rhoref,zz,z)

    ! global dzz_3d kbm2 dz_3d zz_3d;
    real(kind=8), intent(in) :: dti2
    real(kind=8) ::  a1=0.92, b1=16.6, a2=0.74, b2=10.1, c1=0.08, const1
    real(kind=8) ::  e1=1.8, e2=1.33, sef=1.0, cbcnst=100., surfl=2.e5, shiw=0.0
    type(Matrix) ::  l0, kn, boygr, gh,la, dh_3d2, dh, dh_3d,cc,p,utau2
    type(Matrix) :: filter, l0_kappa, tmp, tmp_d
    real(kind=8) :: ghc, coef1, coef2, coef3, coef4, coef5
    type(Matrix) :: LAA, LBB, sh, sm, temp, temp1, temp2
    integer :: ierr
    
    l0    = dm_zeros(im,jm,1);
    kn    = dm_zeros(im,jm,kb);
    boygr = dm_zeros(im,jm,kb);
    gh    = dm_zeros(im,jm,kb);
    !stf   = dm_zeros(im,jm,kb);
    la    = dm_zeros(1, 1, kb);

    dh = h + etf;
    dh_3d=dm_rep(dh,1,1,kb);

    a = dm_zeros(im,jm,kb);
    c = dm_zeros(im,jm,kb);
    !d = dm_zeros(im,jm,kb);


    ! a(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,2:kbm1);
    ! c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2).*dz_3d(:,:,1:kbm2);
    a = SHIFT(dzz_3d, 3, -1) .em. dz_3d
    c = dzz_3d .em. dz_3d

    
    ! a= DIVISION(-dti2 .* (AZF(kq)+umol) , (a .* dh_3d .* dh_3d));
    ! c= DIVISION(-dti2 .* (AZB(kq)+umol) , (c .* dh_3d .* dh_3d));
    dh_3d2 = dm_squ(dh_3d)
    a = (-dti2 * (AZF(kq) + umol)) .ed. (a .em. dh_3d2)
    c = (-dti2 * (AZB(kq) + umol)) .ed. (c .em. dh_3d2)


    ! a(:,:,1)=0.e0;         c(:,:,1)=0.e0;       
    ! a(:,:,kb)=0.e0;        c(:,:,kb)=0.e0;
    a = a .em. REV_MASK_ZZ
    c = c .em. REV_MASK_ZZ

    ! %-----------------------------------------------------------------------
    ! %     The following section solves the equation:
    ! %       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
    ! %     Surface and bottom boundary conditions:
    const1=(16.6e0**(2.e0/3.e0))*sef;

    
    ! % initialize fields that are not calculated on all boundaries
    ! % but are later used there
    ! l0(:,jm)  =0; l0(im,:)=0;
    ! kn(:,:,:)=0;

    ! utau2 = sqrt( AXF(wusurf).^2 +AYF(wvsurf).^2 );
    utau2 = dm_sqrt(dm_squ(AXF(wusurf)) + dm_squ(AYF(wvsurf)));

    
    !  % Surface length scale following Stacey (1999).
    ! l0 = surfl*utau2/grav;
    l0 = (surfl/grav) * utau2


    
    ! uf(:,:,kb) = sqrt( AXF(wubot).^2 +AYF(wvbot).^2 ) .* const1;
    uf = uf .em. REV_MASK_Z2 &
         + dm_rep(const1 * dm_sqrt(dm_squ(AXF(wubot)) + dm_squ(AYF(wvbot))), 1, 1, kb) &
         .em. REV_MASK_Z1

    
    ! %    Calculate speed of sound squared:
    ! h_3d=repmat(h,1,1,kb);
    h_3d = dm_rep(h, 1,1,kb)

    !%     Calculate pressure in units of decibars:
    ! p=grav*rhoref*(-zz_3d .* h_3d)*1.e-4;
    p = -grav*rhoref*1.e-4*(zz_3d .em. h_3d)
    
    ! cc=1449.10+.00821*p+4.55*(t+tbias) -.045e0*(t+tbias).^2 +1.34*(s+sbias-35.0e0);
    cc = 1449.10+0.00821*p+4.55*(t+tbias)-.045e0*dm_squ(t+tbias) + &
         1.34*(s+sbias-35.0e0)
    
    ! cc=cc./sqrt((1.0-.01642.*p./cc) .*(1.0-0.40.*p./cc.^2));
    cc = cc .ed. dm_sqrt((1.0-0.01642*p .ed. cc) .em. (1.0-0.40*p .ed. dm_squ(cc)))
    ! cc(:,:,kb)=0;
    cc = cc .em. REV_MASK_Z2

    
    ! %     Calculate buoyancy gradient:
    ! q2b =abs(q2b);
    q2b = dm_abs(q2b)
    
    ! q2lb=abs(q2lb);
    q2lb = dm_abs(q2lb)
    
    ! tmp = zeros(im,jm,kb);
    ! for k=2:kbm1
    !     tmp(:,:,k)=dzz(k-1);
    ! end
    tmp = SHIFT(dzz_3d, 3, -1)
    
    ! boygr=DIVISION(-grav* DZB(rho) , tmp .* h_3d ) + DIVISION(grav^2 , AZB(cc.^2));
    boygr = (-grav*DZB(rho)) .ed. (tmp .em. h_3d) + grav**2 * dm_pow(AZB(dm_squ(cc)), -1)
    
    ! boygr(:,:,1)=0.e0; boygr(:,:,kb)=0.e0;
    boygr = boygr .em. REV_MASK_ZZ
    
    ! l=q2lb ./ q2b;
    l = q2lb .ed. q2b
    
    ! %l=max(l, repmat(kappa*l0,1,1,kb));
    ! l=(z_3d>-0.5) .* max(l, repmat(kappa*l0,1,1,kb))+(z_3d<=-0.5) .* l;
    l0_kappa = dm_rep(kappa * l0, 1, 1, kb)
    filter = (l > l0_kappa)
    l = (z_3d > -0.5) .em. (filter .em. l + (1 - filter) .em. l0_kappa) + &
         (z_3d <= -0.5) .em. l

    
    ! gh=l.^2 .* boygr ./q2b;
    ! gh=min(gh, 0.028);
    gh = dm_squ(l) .em. boygr .ed. q2b
    filter = (gh >= 0.028)
    gh = 0.028 * filter + gh .em. (1-filter)
    
    ! l(:,:,1)=kappa*l0; l(:,:,kb)=0;
    ! gh(:,:,1)=0      ; gh(:,:,kb)=0;
    l = l .em. REV_MASK_ZZ
    gh = gh .em. REV_MASK_ZZ
    
    ! %    Calculate production of turbulent kinetic energy:
    ! kn= DIVISION(km.*sef.*(AXF(DZB(u)).^2 + AYF(DZB(v)).^2) , ...
    !             (tmp.*dh_3d).^2) -shiw.*km.*boygr + kh.*boygr;
    km = (sef * km .em. (dm_squ(AXF(DZB(u))) + dm_squ(AYF(DZB(v))))) .ed. &
         (dm_squ(tmp .em. dh_3d) - shiw * km .em. boygr + kh .em. boygr)


    ! %  NOTE: Richardson # dep. dissipation correction (Mellor: 2001; Ezer, 2000),
    ! %  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
    ghc=-6.0e0;
    !stf=dm_ones(im,jm,kb);
    
    ! % It is unclear yet if diss. corr. is needed when surf. waves are included.
    ! %           if(gh(i,j,k).lt.0.e0)
    ! %    ...        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
    ! %           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
    ! dtef=sqrt(q2b).*stf./(b1.*l+small);
    ! dtef(:,:,1)=0.e0;dtef(:,:,kb)=0.e0;
    dtef = dm_sqrt(q2b) .ed. (b1 * l + small) .em. REV_MASK_ZZ

    !     d=-uf-2.e0*dti2*kn;
    !     d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
    !     d(:,jm,1)=0; d(im,:,1)=0; d(:,:,kb)=-uf(:,:,kb);
    !     temp1=a(:,:,1:kbm1);    temp2=c(:,:,2:kb);
    tmp_d = -2.0*dti2*kn-uf
    tmp_d = tmp_d .em. REV_MASK_Z1 - (15.8*cbcnst)**(2./3.) * dm_rep(utau2, 1, 1, kb)
    tmp_d = (1-(MASK_Y2 .em. MASK_Z1)) .em. tmp_d
    tmp_d = (1-(MASK_X2 .em. MASK_Z1)) .em. tmp_d
    tmp_d = (REV_MASK_Z1 .em. tmp_d) - (uf .em. MASK_Z1)
    ! call dm_finalize(ierr)
    ! stop
    
    !   for j=2:jm
    !       for i=2:im
    !    la=diag(reshape(a(i,j,:)+c(i,j,:)-1-2.e0*dti2.*dtef(i,j,:),kb,1),0) ...
    !       - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
    !       - diag(reshape(temp2(i,j,:),kbm1,1),-1);
    !    uf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
    !       end
    !   end
    
    !(m,n,k) => (k, k, m*n) 
    LAA = dm_trid(temp2, a+c-1-2.0*dti2*dtef, temp1)
    !(m,n,k) => (k, 1, m*n)
    LBB = dm_trid1(tmp_d)
    uf = dm_trid2(dm_solve(LAA, LBB), im, jm, kb)
    
    ! %-----------------------------------------------------------------------
    ! %     The following section solves the equation:
    ! %
    ! %       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
    ! vf(:,:,kb)=0.e0;
    vf = vf .em. MASK_Z2

    ! for k=2:kbm1
    !     dtef(:,:,k)=dtef(:,:,k).*(1.e0+e2.*((1.e0/abs(z(k)-z(1))+1.e0./abs(z(k)-z(kb)))     ...
    !                         .*l(:,:,k)./(dh(:,:)*kappa)).^2);
    ! end
    !     d=-vf-dti2*kn.*l.*e1;
    !     d(:,:,1)=-(15.8*cbcnst)^(2./3.).*utau2;
    !     d(:,jm,1)=0.e0; d(im,:,1)=0.e0; d(:,:,kb)=0.e0;
    !     temp=vf(:,:,1);
    tmp_d = -dti2 * kn .em. l .em. e1 - vf
    tmp_d = -(15.8*cbcnst)**(2./3.) * utau2
    tmp_d = tmp_d .em. REV_MASK_X1
    tmp_d = tmp_d .em. REV_MASK_Y1
    tmp_d = tmp_d .em. REV_MASK_Z1
    temp = vf .em. MASK_Z1
    
    !   for j=2:jm
    !       for i=2:im
    !    la=diag(reshape(a(i,j,:)+c(i,j,:)-1-dti2.*dtef(i,j,:),kb,1),0) ...
    !       - diag(reshape(temp1(i,j,:),kbm1,1),1) ...
    !       - diag(reshape(temp2(i,j,:),kbm1,1),-1);
    !    vf(i,j,:)=la\reshape(d(i,j,:),kb,1); 
    !       end
    !   end
    LAA = dm_trid(temp2, a+c-1-dti2 * dtef, temp1)
    LBB = dm_trid1(tmp_d)
    vf = dm_trid2(dm_solve(LAA, LBB), im, jm, kb)

    !    vf(:,:,1)=temp;
    vf = vf .em. REV_MASK_Z1 + temp
    
    ! dt_3d=repmat(dt,1,1,kb);
    dt_3d = dm_rep(dt, 1, 1, kb)
    
    ! filter = (uf<=small | vf<=small);
    ! filter(:,:,1) = false;
    ! filter(:,:,kb) = false;
    ! uf(filter) = small;
    ! vf(filter) = 0.1 * dt_3d(filter) * small;
    filter = 1-((1 - (uf <= small)) .em. (1 - (vf <= small)))
    filter = filter .em. REV_MASK_Z1
    filter = filter .em. REV_MASK_Z2
    uf = (uf .em. (1-filter)) + (small * filter)
    vf = 0.1 * small * (dt_3d .em. filter) + (vf .em. (1-filter))
    
    ! %-----------------------------------------------------------------------
    ! %     The following section solves for km and kh:
    !     coef4=18.e0*a1*a1+9.e0*a1*a2;
    !     coef5=9.e0*a1*a2;
    ! %     Note that sm and sh limit to infinity when gh approaches 0.0288:
    !     coef1=a2*(1.e0-6.e0*a1/b1*stf);
    !     coef2=3.e0*a2*b2/stf+18.e0*a1*a2;
    !     coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf);
    !     sh=coef1./(1.e0-coef2.*gh);
    !     sm=coef3+sh.*coef4.*gh;
    !     sm=sm./(1.0-coef5.*gh);
    coef4=18.e0*a1*a1+9.e0*a1*a2;
    coef5=9.e0*a1*a2;
    coef1=a2*(1.e0-6.e0*a1/b1);
    coef2=3.e0*a2*b2+18.e0*a1*a2;
    coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1);
    
    sh=coef1 * dm_pow(1.e0-coef2*gh, -1);
    sm=coef3 + coef4 * sh .em. gh;
    sm=sm .em. dm_pow(1.0-coef5 * gh, -1);

    !     kn=l.*sqrt(abs(q2));
    !     kq=(kn.*.41e0.*sh+kq)*.5e0;
    !     km=(kn.*sm+km)*.5e0;
    !     kh=(kn.*sh+kh)*.5e0;
    kn = l .em. dm_sqrt(dm_abs(q2))
    kq = (0.41*kn .em. sh + kq) * 0.5
    km = (kn .em. sm + km) * 0.5
    kh = (kn .em. sh + kh) * 0.5
    
    !       fsm_3d = repmat(fsm, 1, 1, kb);
    !       km(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
    !       kh(:,jm,:)=km(:,jmm1,:).*fsm_3d(:,jm,:);
    !       km(:,1,:)=km(:,2,:).*fsm_3d(:,1,:);
    !       kh(:,1,:)=kh(:,2,:).*fsm_3d(:,1,:);
    !fsm_3d = dm_rep(fsm, 1, 1, kb)

    !       km(im,:,:)=km(imm1,:,:).*fsm_3d(im,:,:);
    !       kh(im,:,:)=kh(imm1,:,:).*fsm_3d(im,:,:);
    !       km(1,:,:)=km(2,:,:).*fsm_3d(1,:,:);
    !       kh(1,:,:)=kh(2,:,:).*fsm_3d(1,:,:);
    
    km = (km .em. REV_MASK_X2) + (MASK_X2 .em. SHIFT(km, 1, -1) .em. fsm_3d)        
    kh = (kh .em. REV_MASK_X2) + (MASK_X2 .em. SHIFT(kh, 1, -1) .em. fsm_3d)    
    km = km .em. REV_MASK_X1 + MASK_X1 .em. SHIFT(km, 1, 1) .em. fsm_3d    
    kh = kh .em. REV_MASK_X1 + MASK_X1 .em. SHIFT(kh, 1, 1) .em. fsm_3d
    !  end

  end subroutine new_profq

  subroutine new_profu(dti2)
    use dm_op
    use dm
    use grid
    use input
    implicit none
    type(Matrix) :: dh,dh_3d
    real, intent(in) :: dti2
    integer :: i, ierr
    ! function [uf,wubot] = new_profu(uf,etf,h,km,wusurf,cbc,ub,vb)
    ! global im  jm kb dz_3d dzz_3d kbm1 dti2 umol kbm2 dz dum
    
    ! dh = AXB(h+etf);dh(1,:)=1.e0;dh(:,1)=1.e0;
    dh = AXB(h+etf);
    call dm_setvalues(dh, (/0/), (/(i,i=0,jm)/), (/0/), (/(1,i=0,jm)/), ierr)
    call dm_setvalues(dh, (/(i,i=0,im)/), (/0/), (/0/), (/(1,i=0,im)/), ierr)
    
    ! dh_3d=repmat(dh,1,1,kb);
    dh_3d = dm_rep(dh, 1, 1, kb)
    
    ! a=zeros(im,jm,kb);
    ! c = AXB(km); c(1,:,:)=0.e0; c(:,1,:)=0.e0;
    c = AXB(km) .em. REV_MASK_X1 .em. REV_MASK_Y1;

    ! la=zeros(kbm1);d=zeros(im,jm,kb);
    
    ! d(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
    ! a=DIVISION(d,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
    ! a(:,:,kbm1)= 0.e0;
    tmp_d = SHIFT(-dti2 * (c + umol), 3, 1)
    a = tmp_d .ed. (dz_3d .em. dzz_3d .em. dh_3d .em. dh_3d)
    a = (1 - SHIFT(MASK_Z2, 3, 1)) .em. a
    
    ! d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    ! c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);
    ! c(:,:,1)  =0.e0;
    tmp_d = (REV_MASK_ZZ .em. tmp_d) + (SHIFT(dzz_3d, 3, -1) .em. REV_MASK_Z2)
    c = (-dti2*(c+umol)) .ed. (dz_3d .em. tmp_d .em. dh_3d .em. dh_3d)
    c = c .em. REV_MASK_Z1
    
    ! tps = AXB(cbc) .* sqrt( ub(:,:,kbm1).^2 + AXB( AYF( vb(:,:,kbm1) ) ).^2 );
    
    !     d=-uf;
    !     d(:,:,1)= -uf(:,:,1) + dti2 .* wusurf(:,:) ./ (dh(:,:) .* dz(1));
    !     d(:,:,kbm1)=-uf(:,:,kbm1)+tps.*dti2./(dz(kbm1)*dh); 
    
    !   for j=2:jm
    !       for i=2:im
    !    la=diag(reshape(a(i,j,1:kbm1)+c(i,j,1:kbm1)-1,kbm1,1),0) ...
    !       - diag(reshape(a(i,j,1:kbm2),kbm2,1),1) ...
    !       - diag(reshape(c(i,j,2:kbm1),kbm2,1),-1);
    !    uf(i,j,1:kbm1)=la\reshape(d(i,j,1:kbm1),kbm1,1); 
    !       end
    !   end
    !    uf=uf.*repmat(dum,1,1,kb);
    !    wubot=-tps.*uf(:,:,kbm1);
    
    ! return
  end subroutine

  ! subroutine new_profv()
    
  ! end subroutine new_profu  
end module functions
