program pom2k
  use dm_op 
  use dm
  use dm_type
  use input
  use grid
  
  implicit none
  type(Matrix) :: tmp
  integer :: ierr
  integer :: myrank, mysize
  real(kind=8) :: tmp_val, val1, val2
  integer      :: pos1(3), pos2(3)
  real(kind=8) :: array_one(1)
  real(kind=8) :: time0, time, period
  real(kind=8) :: dti, dte2, dti2, ispi, isp2i
  integer, allocatable :: idxm1(:), idxn1(:), array1(:)

  type(Matrix) :: tmp_dvol, darea
  type(Matrix) :: dt_3d1, dt_axb, dt_ayb, tps_3d, utb_3d, utf_3d, vtf_3d, vtb_3d
  
  real(kind=8) :: vtot, atot, taver, saver, eaver, tsalt
  integer :: i, iprint, iswtch, iexit, iext, iint, imax, jmax
  real(kind=8) :: vamax

  call dm_init1(ierr)
  call dm_comm_rank(myrank,ierr)
  call dm_comm_size(mysize,ierr)

  call InitGrid()
  call InitOperatorModule(im, jm, kb)
  
  call dm_init2(im, jm , kb,ierr)

  allocate(idxm1(imm2), idxn1(jmm2))
  idxm1=(/(i, i=1,imm2)/)
  idxn1=(/(i, i=1,jmm2)/)
  array1=(/(0,i=1,imm2*jmm2)/)

  call file2ic()
  
  ! call dm_view(z, ierr)
  ! call dm_view(zz, ierr)
  ! call dm_view(dzz, ierr)
  ! call dm_view(dz, ierr)  
  
  if(myrank==0) then
     print *, "==============Input paramenters==========="
     print *, "im=",im,",jm=",jm,",kb=",kb
  endif

  !*************************************
  dti=dte * isplit;

  dte2 = dte * 2;
  dti2 = dti * 2;
  iend=max(floor(days*24.0*3600.0/dti+0.5),2);

  iprint=floor(prtd1*24.0*3600.0/dti+0.5);
  iswtch=floor(swtch*24.0*3600.0/dti+0.5);

  ispi=1.0/isplit;
  isp2i=1.0/(2.0*isplit);

  if(myrank==0) then
     print*, "dti=",dti,"dte2=",dte2,"dti2=",dti2,"iend=",iend
     print*, "iprint=",iprint,"iswtch=",iswtch
     print*, "ispi=",ispi,"isp2i=",isp2i
  endif

  !call new_depth(z, zz, dz, dzz, z_3d, zz_3d, dz_3d, dzz_3d, kl1, kl2, ierr)
  call depth(ierr)
  
  !********************
  !LOAD GRIDFILE2IC HERE
  !*******************
  !Inertial period for temporal filter:
  !period=(2.0*pi)/abs(cor(floor(im/2),floor(jm/2)))/86400.0;
  tmp_val = dm_getvalue(cor, im/2, jm/2, 0)
  if(tmp_val == 0.) tmp_val = 1.0
  period=(2.0*pi)/abs(tmp_val)/86400.0;

  !%     Initialise time:
  time0 =0.0;
  time = 0.0;

  ! %     Initial conditions:
  ! %     NOTE that lateral thermodynamic boundary conditions are often set
  ! %     equal to the initial conditions and are held constant thereafter.
  ! %    Users can of course create variable boundary conditions.
  ua = uab;
  va = vab;
  el = elb;
  et = etb;
  etf= et;
  d = h + el;
  dt = h + et;

  !w(:,:,1)=vfluxf;
  w = w - w .em. MASK_Z1 + dm_rep(vfluxf, 1, 1, kb) .em. MASK_Z1

  ! d_3d=repmat(d,1,1,kb);
  ! dt_3d=repmat(dt,1,1,kb);
  d_3d  = dm_rep(d,1,1,kb);
  dt_3d = dm_rep(dt,1,1,kb);

  l    = 0.1 * dt_3d;
  q2b  = small * ONES;
  q2lb        = l .em. q2b;
  kh          = l .em. dm_sqrt(q2b);
  km          = kh;
  kq          = kh;
  aam = aam_init * ONES;

  q2 = q2b;
  q2l= q2lb;
  t  = tb;
  s  = sb;
  u  = ub;
  v  = vb;
  
  call dense()

  call baropg()  
  ! drx2d = sum(drhox .* dz_3d .* REV_MASK_KB, 3);
  ! dry2d = sum(drhoy .* dz_3d .* REV_MASK_KB, 3);
  drx2d = CSUM(drhox .em. dz_3d .em. REV_MASK_Z2, 4)
  dry2d = CSUM(drhoy .em. dz_3d .em. REV_MASK_Z2, 4)

  ! % Calculate bottom friction coefficient:
  ! cbc=(kappa./log((1.0+zz(kbm1))*h/z0b)).^2;
  ! cbc=max(cbcmin,cbc);

  cbc=dm_pow(kappa * dm_pow(dm_log((1.0+dm_getvalue(zz, 0, 0, kbm1-1)) * (1.0/z0b) * h), -1.0), 2.0);
  !!cbc=max(cbcmin,cbc);

  ! % If the following is invoked, then it is probable that the wrong
  ! % choice of z0b or vertical spacing has been made:

  !! cbc=min(cbcmax,cbc);

  ! % Calculate external (2-D) CFL time step:
  ! tps=0.5./sqrt(1.0./dx.^2+1.0./dy.^2) ./ sqrt(grav*(h+small)) .* fsm;
  tps=0.5 * dm_pow((dm_sqrt(dm_pow(dx,-2)+dm_pow(dy,-2)) .em. &
       dm_sqrt(grav*(h+small))), -1) .em. fsm;

  d = h + el;
  dt = h + et;
  time = time0;

  d_3d=dm_rep(d,1,1,kb);
  dt_3d=dm_rep(dt,1,1,kb);

  ! call dm_finalize(ierr)
  ! stop
  ! %==========================================
  ! %           begin internal (3-D) mode
  ! %==========================================
  !for iint=1:iend
  do iint = 1, iend
     time=dti*iint*1.0/86400.0+time0;
     if(lramp) then
        ramp = time/period;
        if(ramp>1.0) ramp=1.0;
     else
        ramp=1.0;
     endif
     
     ! %       write(6,2) mode,iint,time
     ! %   2   format(' mode,iint,time =',2i5,f9.2)
     ! %-----------------------------------------------------------------------
     ! %     Set time dependent, surface and lateral boundary conditions.
     ! %     The latter will be used in subroutine bcond. Users may
     ! %     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
     ! %     wssurf, swrad and vflux.
     ! %
     ! %     Introduce simple wind stress. Value is negative for westerly or
     ! %     southerly winds. The following wind stress has been tapered
     ! %     along the boundary to suppress numerically induced oscilations
     ! %     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
     ! %     To make a healthy surface Ekman layer, it would be well to set
     ! %     kl1=9.
     ! %
     call dm_setvalues(e_atmos, idxm1, idxn1, (/0/), array1, ierr)
     call dm_setvalues(vfluxf, idxm1, idxn1, (/0/), array1, ierr)
     call dm_setvalues(w, idxm1, idxn1, (/0/), array1, ierr)
     call dm_setvalues(swrad,  idxm1, idxn1, (/0/), array1, ierr)    
     call dm_setvalues(wtsurf, idxm1, idxn1, (/0/), array1, ierr)
     call dm_setvalues(wssurf, idxm1, idxn1, (/0/), array1, ierr)

     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__

     satm = 0.0
  
     ! %-----------------------------------------------------------------------
     ! %
     ! %     Set lateral viscosity:
     ! %
     ! %     If mode=2 then initial values of aam2d are used. If one wishes
     ! %     to use Smagorinsky lateral viscosity and diffusion for an
     ! %     external (2-D) mode calculation, then appropiate code can be
     ! %     adapted from that below and installed just before the end of the
     ! %     "if(mode.eq.2)" loop in subroutine advave.
     ! %
     ! %     %alculate Smagorinsky lateral viscosity:
     ! %
     ! %       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
     ! %                                     +.5*(du/dy+dv/dx)**2) )
     ! %
      
     if(mode/=2) then
        ![advx,advy]=new_advct(u,v,dt_3d,aam,ub,vb);
        call advct()
     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        ! [rho,drhox,drhoy] = new_baropg(rho, rmean, dt, ramp);
        call baropg()
     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        ! aam=horcon .* dx_3d .* dy_3d .*sqrt((DXF(u)./dx_3d).^2 + (DYF(v)./dy_3d).^2 ...
        ! +0.5*( DYB(AYF(AXF(u)))./dy_3d + DXB(AXF(AYF(v)))./dx_3d).^2 );
        aam = horcon * dx_3d .em. dy_3d .em. dm_sqrt(dm_squ(DXF(u) .ed. dx_3d) + &
             dm_squ((DYF(v) .ed. dy_3d)) + 0.5 * dm_squ(DYB(AYF(AXF(u))) .ed. dy_3d + &
             DXB(AXF(AYF(v))) .ed. dx_3d))

     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        ! aam(:,1,:)=aam_init;
        ! aam(:,jm,:)=aam_init;
        ! aam(1,:,:)=aam_init;
        ! aam(im,:,:)=aam_init;
        ! aam(:,:,kb)=aam_init;
        aam = aam - (aam .em. MASK_Y1 - aam_init * MASK_Y1)
        aam = aam - (aam .em. MASK_Y2 - aam_init * MASK_Y2)
        aam = aam - (aam .em. MASK_X1 - aam_init * MASK_X1)
        aam = aam - (aam .em. MASK_X2 - aam_init * MASK_X2)
        aam = aam - (aam .em. MASK_Z2 - aam_init * MASK_Z2)

     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        ! % Form vertical averages of 3-D fields for use in external (2-D)
        ! % mode:
        ! adx2d = sum(advx.*dz_3d, 3);
        ! ady2d = sum(advy.*dz_3d, 3);
        ! drx2d = sum(drhox.*dz_3d, 3);
        ! dry2d = sum(drhoy.*dz_3d, 3);
        ! aam2d = sum(aam.*dz_3d, 3);
        adx2d = CSUM(advx .em. dz_3d, 4)
        ady2d = CSUM(advy .em. dz_3d, 4);
        drx2d = CSUM(drhox .em. dz_3d, 4);
        dry2d = CSUM(drhoy .em. dz_3d, 4);
        aam2d = CSUM(aam .em. dz_3d, 4);

        if (myrank == 0) &
             print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__

        ! [tps,wubot,wvbot,advua,advva]
        !      = new_advave(tps,wubot,wvbot,mode,aam2d,uab,vab,ua,va,cbc,d);
        call advave(iint)
        
     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__

        adx2d = adx2d - advua;
        ady2d = ady2d - advva;
     endif

     ! call dm_finalize(ierr)
     ! stop
     
     egf= ispi * el;
     
     ! utf=ua .* 2.0 .* AXB(d) .* isp2i;
     ! vtf=va .* 2.0 .* AYB(d) .* isp2i;
     utf = 2.0 * isp2i * ua .em. AXB(d)
     vtf = 2.0 * isp2i * va .em. AYB(d)
     
     !% Begin external (2-D) mode        
     do iext=1,isplit    
        ! %     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
        ! %     with pom98.f. See also modifications to subroutine vertvl.
        !elf= elb+dte2.*(-(DXF( AXB(d).*AXB(dy).*ua)+DYF(AYB(d).*AYB(dx).*va))./ art-vfluxf);  
        elf= elb + dte2 * (-(DXF(AXB(d) .em. AXB(dy) .em. ua) + &
             DYF(AYB(d) .em. AYB(dx) .em. va)) .em. dm_pow(art, -1) - vfluxf)

        ! [elf,uaf,vaf,uf,vf,w] = new_bcond(1,elf,uaf,vaf,uf,vf,w,...
        !     im,jm,kb,imm1,jmm1,kbm1,...
        !     fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
        !     dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        !call bcond(1, dti)
        if(mod(iext,ispadv)==0) then
           ![tps,wubot,wvbot,advua,advva] = new_advave(tps,wubot,wvbot,mode,aam2d,uab,vab,ua,va,cbc,d);
           call advave(iint)
        endif

     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        ! uaf=   DIVISION( (2.0*AXB(h+elb) .* aru .* uab ...
        !        -4.0* dte .* (adx2d + advua - aru .* AXB(cor .* d .* AYF(va)) ...
        !        + grav .* AXB(dy) .* AXB(d)  ...
        !          .*( (1.0-2.0*alpha) .* DXB(el) + alpha* (DXB(elb)+ DXB(elf)) + DXB(e_atmos) ) ...
        !        + drx2d + aru .* (wusurf-wubot))) , (2.0*AXB(h+elf) .* aru));      
        uaf = (2.0*AXB(h+elb) .em. aru .em. uab &
             -4.0 * dte .em. (adx2d + advua - aru .em. AXB(cor .em. d .em. AYF(va)) + &
             grav * AXB(dy) .em. AXB(d) .em. &
             ((1.0-2.0*alpha) * DXB(el) + alpha * (DXB(elb)+ DXB(elf)) + DXB(e_atmos) ) + &
             drx2d + aru .em. (wusurf-wubot))) .ed. (2.0*AXB(h+elf) .em. aru);      

        vaf = (2.0*AYB(h+elb) .em. arv .em. vab &
             -4.0* dte .em. (ady2d + advva + arv .em. AYB(cor .em. d .em. AXF(ua)) + &
             grav * AYB(dx) .em. AYB(d) .em. &
             ((1.0-2.0*alpha) * DYB(el) + alpha * (DYB(elb)+ DYB(elf)) + DYB(e_atmos)) + &
             dry2d + arv .em. (wvsurf-wvbot))) .ed. (2.0*AYB(h+elf) .em. arv);  

        ! [elf,uaf,vaf,uf,vf,w] = new_bcond(2,elf,uaf,vaf,uf,vf,w,...
        !     im,jm,kb,imm1,jmm1,kbm1,...
        !     fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
        !     dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        !call bcond(2, dti)
        ! if(iext==(isplit-2))
        !     etf=0.25*smoth*elf;
        ! elseif(iext==(isplit-1))
        !     etf=etf+0.5*(1.0-0.5*smoth)*elf;
        ! elseif(iext==isplit)
        !     etf=(etf+0.5*elf).*fsm;
        ! end
        if(iext==(isplit-2)) then
           etf=0.25*smoth*elf;
        elseif(iext==(isplit-1)) then
           etf=etf+0.5*(1.0-0.5*smoth)*elf;
        elseif(iext==isplit) then
           etf=(etf+0.5*elf) .em. fsm;
        endif

     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
      
        ! % Stop if velocity condition violated (generally due to CFL
        ! % criterion not being satisfied):   
        ! [vamax, vapos]=max(abs(vaf(:)));
        ! [imax, jmax]=ind2sub(size(vaf), vapos);
        call dm_max(vaf, val1, pos1, ierr)
        call dm_min(vaf, val2, pos2, ierr)

        vamax = max(abs(val1), abs(val2))

        ! if(vamax <= vmaxl)
        !     ! %
        !     ! %     Apply filter to remove time split and reset time sequence:
        !     ! %
        !     ua=ua+0.5*smoth*(uab-2.0*ua+uaf);
        !     va=va+0.5*smoth*(vab-2.0*va+vaf);
        !     el=el+0.5*smoth*(elb-2.0*el+elf);
        !     elb=el;
        !     el=elf;
        !     d=h+el;
        !     uab=ua;
        !     ua=uaf;
        !     vab=va;
        !     va=vaf;
        !     if(iext~=isplit)
        !         egf=egf+el*ispi;
        !         utf=utf+2.0* ua .* AXB(d) * isp2i;
        !         vtf=vtf+2.0* va .* AYB(d) * isp2i;
        !     end
        !  end do

        if(vamax <= vmaxl) then
           ua = ua+0.5*smoth*(uab-2.0*ua+uaf);
           va = va+0.5*smoth*(vab-2.0*va+vaf);
           el = el+0.5*smoth*(elb-2.0*el+elf);
           elb = el;
           el = elf;
           d = h+el;
           uab = ua;
           ua = uaf;
           vab = va;
           va = vaf;

           if(iext /= isplit) then
              egf = egf + el * ispi
              utf = utf + 2.0 * ua * isp2i .em. AXB(d)
              vtf = vtf + 2.0 * va * isp2i .em. AYB(d)
           endif
        endif
     enddo

     if (myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
     
     ! call dm_finalize(ierr)
     ! stop
     ! %===========================================
     ! %End of external (2-D) mode
     ! %=============================================
     if(vamax <= vmaxl) then
        ! continue with internal (3-D) mode calculation:
        if((iint /= 1 .or. time0 /= 0.e0) .and. mode/=2) then

           !% Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
           !tps=sum(u(:,:,1:kbm1).*dz_3d(:,:,1:kbm1), 3);
           tps=CSUM(u .em. dz_3d .em. MASK_Z2, 4)

           ! utb_3d = repmat(utb, 1, 1, kb);
           ! utf_3d = repmat(utf, 1, 1, kb);
           ! tps_3d = repmat(tps, 1, 1, kb);
           ! dt_3d1 = repmat(dt,  1, 1, kb);
           utb_3d = dm_rep(utb, 1, 1, kb);
           utf_3d = dm_rep(utf, 1, 1, kb);
           tps_3d = dm_rep(tps, 1, 1, kb);
           dt_3d1 = dm_rep(dt,  1, 1, kb);

           !dt_axb = 2.0 * AXB(dt_3d1);
           dt_axb = 2.0 * AXB(dt_3d1);

           !u = u - (tps_3d - (utb_3d+utf_3d) ./ dt_axb) .* REV_MASK_KB;
           u = u - (tps_3d - (utb_3d+utf_3d) .ed. dt_axb) .em. REV_MASK_Z2;

           !tps = sum(v(:,:,1:kbm1) .* dz_3d(:,:,1:kbm1), 3);
           tps = CSUM(v .em. dz_3d .em. MASK_Z2, 4);

           ! vtb_3d = repmat(vtb, 1, 1, kb);
           ! vtf_3d = repmat(vtf, 1, 1, kb);
           ! tps_3d = repmat(tps, 1, 1, kb);
           vtb_3d = dm_rep(vtb, 1, 1, kb);
           vtf_3d = dm_rep(vtf, 1, 1, kb);
           tps_3d = dm_rep(tps, 1, 1, kb);

           !dt_ayb = 2.0 * AYB(dt_3d1);                      
           dt_ayb = 2.0 * AYB(dt_3d1);

           !v = v - (tps_3d - (vtb_3d+vtf_3d) ./ dt_ayb) .* REV_MASK_KB;           
           v = v - (tps_3d - (vtb_3d+vtf_3d) .ed. dt_ayb) .em. REV_MASK_Z2;

           ![a0,c0,w0]=new_vertvl(a,c,w,dx,dy,dz,dt,u,v,vfluxb, ...
           !          vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1);
           call vertvl(dti2)
           
           ! [elf,uaf,vaf,uf,vf,w] = new_bcond(5,elf,uaf,vaf,uf,vf,w,...
           !  im,jm,kb,imm1,jmm1,kbm1,...
           !  fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
           !  dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs, ...
           !  sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
           !call bcond(5, dti)
           
           vf = ZEROS
           uf = ZEROS


           ! THIS FUNCTION SHOULD BE REPLACED
           !uf=new_advq(q2b,q2,dt,u,v,w,aam,h,etb,etf,dti2);    
           call advq(uf, q2b, q2, dti2)
           
           ! THIS FUNCTION SHOULD BE REPLACED            
           !vf=new_advq(q2lb,q2l,dt,u,v,w,aam,h,etb,etf,dti2);  
           call advq(vf, q2lb, q2l, dti2)

           if (myrank == 0) &
                print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
           ! call dm_finalize(ierr)
           ! stop
           
           ! THIS FUNCTION SHOULD BE REPLACED
           ! [a,c,tps,dtef,...
           !     ee,gg,l,kq,km,kh,...
           !     uf,vf,q2b,q2lb,a,c]=new_profq(a,c,tps,dtef,....
           !     ee,gg,l,kq,km,kh,...
           !     uf,vf,q2,q2b,q2lb,a,c,...
           !     h,etf,dti2,umol,dzz,grav,rho,kappa,
           !     u,v,dt,small,fsm,im,jm,kb,imm1,jmm1,kbm1,tbias,sbias,dz,...
           !     wusurf,wubot,wvsurf,wvbot,t,s,rhoref,zz,z);
           !call profq(dti2)

           if (myrank == 0) &
                print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
           ! call dm_finalize(ierr)
           ! stop
           
           ! call dm_finalize(ierr)
           ! stop
           ! THIS FUNCTION SHOULD BE REPLACED           
           ! [elf,uaf,vaf,uf,vf,w] = new_bcond(6,elf,uaf,vaf,uf,vf,w,...
           !   im,jm,kb,imm1,jmm1,kbm1,...
           !   fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
           !   dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs, ...
           !   sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
           !call bcond(6, dti)
           
           q2 = q2  + .5e0 * smoth * (uf + q2b - 2.e0 * q2);
           q2l= q2l + .5e0 * smoth * (vf + q2lb - 2.e0 * q2l);
           q2b = q2;
           q2  = uf;
           q2lb = q2l;
           q2l  = vf;

            if (myrank == 0) &
                 print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
           ! call dm_finalize(ierr)
           ! stop
           
           ! calculate tf and sf using uf, vf, a and c as temporary variables:
           if(mode /= 4) then
              if(nadv == 1)  then
                 ! THIS FUNCTION SHOULD BE REPLACED                             
                 !uf=new_advt1(tb,t,dt,u,v,aam,tprni,w,etb,etf,h);
                 !uf = advt1(tb, t, dti2)
                 call advt1(uf, tb, t, dti2)
                 ! THIS FUNCTION SHOULD BE REPLACED                             
                 !vf=new_advt1(sb,s,dt,u,v,aam,tprni,w,etb,etf,h);
                 !vf = advt1(sb, s, dti2)
                 call advt1(vf, sb, s, dti2)
              elseif(nadv==2) then
                 ! THIS FUNCTION SHOULD BE REPLACED                
                 ! [tb,t,tclim,uf,a,c,nitera,sw,...
                 !  zflux] = advt2(tb,t,tclim,uf,a,c,nitera,sw,...
                 !  zflux,...
                 !  im,jm,kb,imm1,jmm1,kbm1,dti2,...
                 !  etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);

                 ! THIS FUNCTION SHOULD BE REPLACED
                 ! [sb,s,sclim,vf,a,c,nitera,sw,...
                 !  zflux] = advt2(sb,s,sclim,vf,a,c,nitera,sw,...
                 !  zflux,...
                 !  im,jm,kb,imm1,jmm1,kbm1,dti2,...
                 !  etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);  
              else
                 print*, "Invalid value for nadv ..... "
                 print*, "program terminated"
                 stop;
              end if

              ! THIS FUNCTION SHOULD BE REPLACED         
              ! [uf] = new_proft(uf,wtsurf,tsurf,nbct,h,etf,swrad,kh);
              ! THIS FUNCTION SHOULD BE REPLACED
              ! [vf] = new_proft(vf,wssurf,ssurf,nbcs,h,etf,swrad,kh);
              call proft(uf, wtsurf, dti2)
              call proft(vf, wssurf, dti2)

            if (myrank == 0) &
                 print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__              
              ! [elf,uaf,vaf,uf,vf,w] = new_bcond(4,elf,uaf,vaf,uf,vf,w,...
              !  im,jm,kb,imm1,jmm1,kbm1,...
              !  fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
              !  dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,
              !  sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
              !call bcond(4, dti)
              
              t = t+.5e0*smoth*(uf+tb-2.e0*t);
              s = s+.5e0*smoth*(vf+sb-2.e0*s);
              tb = t;
              t = uf;
              sb = s;
              s = vf;

              call dense()
               if (myrank == 0) &
                    print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__              
           endif
           
           !% calculate uf and vf:
           ! THIS FUNCTION SHOULD BE REPLACED
           !uf = new_advu(advx,cor,dt,e_atmos,drhox,h,ub,u,v,w,egf,egb,etf,etb);
           ! THIS FUNCTION SHOULD BE REPLACED
           !vf = new_advv(advy,cor,dt,e_atmos,drhoy,h,vb,u,v,w,egf,egb,etf,etb);
           
           call advu(dti2)
           call advv(dti2)
           !vf0 = vf;

            if (myrank == 0) &
                 print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
           
           ! THIS FUNCTION SHOULD BE REPLACED                                            
           ![uf,wubot] = new_profu(uf,etf,h,km,wusurf,cbc,ub,vb);           
           ![vf,wvbot] = new_profv(vf,etf,h,km,wvsurf,cbc,ub,vb);
           call profu(dti2)
           call profv(dti2)

            if (myrank == 0) &
                 print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
           
           ! [elf,uaf,vaf,uf,vf,w] = new_bcond(3,elf,uaf,vaf,uf,vf,w,...
           !    im,jm,kb,imm1,jmm1,kbm1,...
           !    fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
           !    dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,
           !    small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
           !call bcond(3, dti)
           
           !tps = sum((uf+ub-2.e0*u).*dz_3d .* REV_MASK_KB, 3);
           tps = CSUM((uf+ub-2.e0*u) .em. dz_3d .em. REV_MASK_Z2, 4);

           !tps_3d = repmat(tps, 1, 1, kb);
           tps_3d = dm_rep(tps, 1, 1, kb);

           !u = u + 0.5e0*smoth*(uf + ub - 2.e0*u - tps_3d) .* REV_MASK_KB;
           u = u + 0.5e0*smoth*(uf + ub - 2.e0*u - tps_3d) .em. REV_MASK_Z2;

           !tps = sum((vf+vb-2.e0*v) .* dz_3d .* REV_MASK_KB, 3);
           tps = CSUM((vf+vb-2.e0*v) .em. dz_3d .em. REV_MASK_Z2, 4);

           !tps_3d = repmat(tps,1,1,kb);
           tps_3d = dm_rep(tps,1,1,kb);

           !v = v + 0.5e0*smoth*(vf + vb - 2.e0*v-tps_3d) .* REV_MASK_KB;
           v = v + 0.5e0*smoth*(vf + vb - 2.e0*v-tps_3d) .em. REV_MASK_Z2;

           ub = u;
           u = uf;
           vb = v;
           v = vf;
           
            if (myrank == 0) &
                 print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__           
        end if

        egb = egf;
        etb = et;
        et = etf;
        dt = h+et;
        utb = utf;
        vtb = vtf;
        vfluxb = vfluxf;
     end if   !% end if

     ! % Beginning of print section:
     if(iint>=iswtch) then
        iprint=floor(prtd2*24.e0*3600.e0/dti + 0.5);
     endif

      if (myrank == 0) &
           print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__

      if (myrank == 0) &
           print *, "rank=", myrank, " vamax= ", vamax, "vmaxl=", vmaxl
     
     if(mod(iint,iprint)==0 .or. vamax >= vmaxl) then

        if (myrank == 0) &
              print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        if(myrank == 0) then
           print*, "**************************************************"
           write(*, "(A,F15.3,A,I10, A, I10,A, I10)") &
                "time = ", time, "iint = ", iint, &
                "iexit=", iexit, "iprint = ", iprint 
        endif
        
        vtot=0.e0;
        atot=0.e0;
        taver=0.e0;
        saver=0.e0;
        eaver=0.e0;
        
        !dt_3d1 = repmat(dt, 1, 1, kb);
        dt_3d1 = dm_rep(dt, 1, 1, kb)
        
        !tmp_dvol = dx_3d(:,:,1:kbm1).*dy_3d(:,:,1:kbm1).*fsm_3d(:,:,1:kbm1).*
        !           dt_3d1(:,:,1:kbm1).*dz_3d(:,:,1:kbm1);
        tmp_dvol = (dx_3d .em. MASK_Z1) .em. dy_3d .em. fsm_3d .em. dt_3d1 .em. dz_3d;
        
        ! vtot=sum(tmp_dvol(:));
        ! taver=sum(reshape(tb(:,:,1:kbm1).*tmp_dvol,[],1));
        ! saver=sum(reshape(sb(:,:,1:kbm1).*tmp_dvol,[],1));
        vtot  = dm_sum_all(tmp_dvol)
        taver = dm_sum_all(tb .em. tmp_dvol)
        saver = dm_sum_all(sb .em. tmp_dvol)

         if (myrank == 0) &
              print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__
        
        !darea = dx.*dy.*fsm;
        darea = dx .em. dy .em. fsm;

        !atot = sum(darea(:));
        atot = dm_sum_all(darea);

        !eaver = sum(reshape(et.*darea,[],1));
        eaver = dm_sum_all(et .em. darea)

        if(vtot == 0) vtot = 1
        if(atot == 0) atot = 1
        
        taver=taver/vtot;
        saver=saver/vtot;
        eaver=eaver/atot;
        tsalt=(saver+sbias)*vtot;
        
        ! fprintf('vtot = %.6f,atot = %.6f\n',vtot,atot);
        ! fprintf('eaver = %.6f,taver = %.6f,saver=%.6f,saver = %.6f,tsalt = %.6f\n', &
        !      eaver,taver,saver,tsalt,tsalt);

        if(myrank == 0) then
           write(*,*) "vtot = ", vtot, " atot = ", atot 
           write(*,*) "eaver = ", eaver, " taver = ", taver, &
                " saver = ", saver, "tsalt = ",tsalt;
        endif
        
        if(vamax>vmaxl) then
           if(myrank == 0) then
              ! fprintf('time = %.6f, iint = %.6f, iexit = %.6f, iprint = %.6f\n', &
              !      time,iint,iext,iprint);
              write(*,*) "time = ", time, "iint = ", iint, &
                   "iexit = ", iexit, "iprint=", iprint

              ! printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s, &
              !      rho,aam,km,kb,mode,dt,zz,z);

              write(*,*) "************************************************"
              write(*,*) "************ abnormal job end ******************"
              write(*,*) "************* user terminated ******************"
              write(*,*) "************************************************"
              !fprintf('vamax = %d, imax = %d ,jmax = %d \n',vmax,imax,jmax);
              write(*,*) "vamax=",vamax, " imax=", imax, " jmax=", jmax
           endif
           stop
        end if
     end if
     if(myrank == 0) &
          print *, "rank=", myrank, " in file: ", __FILE__, "line:", __LINE__    
     ! %-----------------------------------------------------------------------
     ! %
     ! %  End of internal (3-D) mode
     ! %-----------------------------------------------------------------------  
  end do

  call FinalizeOperatorModule()

  call FinalizeGrid()  
  call dm_finalize(ierr)

end program pom2k
