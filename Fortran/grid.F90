module grid
  use dm_type
  use dm
  use dm_op
  use input
  implicit none

  ! type(Matrix) :: dx, dy, cor, east_c, east_e, east_u, east_v
  ! type(Matrix) :: north_c, north_e, north_u, north_v, h, art, aru
  ! type(Matrix) :: arv, fsm, dum, dvm, tb, sb, tclim, sclim, ub, uab, elb, etb
  type(Matrix) :: dt,    tbe

  real(kind=8)    rfe, rfw, rfn, rfs
  type(Matrix) :: dx_3d, dy_3d, cor_3d, h_3d, art_3d, aru_3d, arv_3d
  type(Matrix) :: fsm_3d, dum_3d, dvm_3d, dt_3d, d_3d


  type(Matrix) aam2d, advua, advva, adx2d, ady2d, art, aru, arv, cbc, cor, d, drx2d
  type(Matrix) dry2d   , dum     , dvm     
  type(Matrix) dx      , dy      , east_c  , east_e  
  type(Matrix) east_u  , east_v  , e_atmos , egb     
  type(Matrix) egf     , el      , elb     , elf     
  type(Matrix) et      , etb     , etf     , fluxua  
  type(Matrix) fluxva  , fsm     , h       , north_c 
  type(Matrix) north_e , north_u , north_v , psi     
  type(Matrix) rot     , ssurf   , swrad   , vfluxb  
  type(Matrix) tps     , tsurf   , ua      , vfluxf  
  type(Matrix) uab     , uaf     , utb     , utf     
  type(Matrix) va      , vab     , vaf    
  type(Matrix) vtb     , vtf     , wssurf  , wtsurf  
  type(Matrix) wubot   , wusurf  , wvbot   , wvsurf  
  type(Matrix) wubot1  , wvbot1  

  type(Matrix) aam  , advx , advy , a
  type(Matrix) c    , drhox, drhoy, dtef  
  type(Matrix) ee   , gg   , kh   , km    
  type(Matrix) kq   , l    , q2b  , q2    
  type(Matrix) q2lb , q2l  , rho  , rmean
  type(Matrix) sb   , sclim, s    , tb    
  type(Matrix) tclim, t    , ub   , uf    
  type(Matrix) u    , vb   , vf   , v    
  type(Matrix) w    

  type(Matrix) ele      , eln      , els   , elw       
  type(Matrix) sbe,  sbn     , sbs     , sbw      
  type(Matrix) tbn     , tbs     , tbw      
  type(Matrix) uabe     , uabw     , ube     , ubw
  type(Matrix) vabn     , vabs     , vbn     , vbs    

  type(Matrix) z, zz, dz, dzz, z_3d, zz_3d, dz_3d, dzz_3d
  type(Matrix) xflux, yflux, zflux, curv, curv2d

  type(Matrix)   :: etb_3d,  egf_3d, egb_3d
  type(Matrix)   :: e_atmos_3d, etf_3d
contains
  subroutine InitGrid()

    call LoadInput()

    xflux = dm_zeros(im,jm,kb);
    yflux = dm_zeros(im,jm,kb);
    curv = dm_zeros(im,jm,kb);
    curv2d = dm_zeros(im, jm, 1)
    
    dx=dm_zeros(im,jm,1)      ;dy=dm_zeros(im,jm,1)      ;
    cor=dm_zeros(im,jm,1)     ;
    east_c=dm_zeros(im,jm,1)  ;east_e=dm_zeros(im,jm,1)  ;
    east_u=dm_zeros(im,jm,1)  ;east_v=dm_zeros(im,jm,1)  ;
    north_c=dm_zeros(im,jm,1) ;north_e=dm_zeros(im,jm,1) ;
    north_u=dm_zeros(im,jm,1) ;north_v=dm_zeros(im,jm,1) ;
    h=dm_zeros(im,jm,1)       ;art=dm_zeros(im,jm,1)     ;
    aru=dm_zeros(im,jm,1)     ;arv=dm_zeros(im,jm,1)     ;
    fsm=dm_zeros(im,jm,1)     ;dum=dm_zeros(im,jm,1)     ;
    dvm=dm_zeros(im,jm,1)     ;
    tb=dm_zeros(im,jm,kb)     ;sb=dm_zeros(im,jm,kb)   ;
    tclim=dm_zeros(im,jm,kb)  ;sclim=dm_zeros(im,jm,kb);
    ub=dm_zeros(im,jm,kb)     ;uab=dm_zeros(im,jm,1)     ;
    elb=dm_zeros(im,jm,1)     ;etb=dm_zeros(im,jm,1)     ;
    dt=dm_zeros(im,jm,1)      ;aam2d=dm_zeros(im,jm,1)   ;
    rho=dm_zeros(im,jm,kb)    ;rmean=dm_zeros(im,jm,kb);
    uabw=dm_zeros(1,jm,1)     ;uabe=dm_zeros(1,jm,1)     ;
    ele=dm_zeros(1,jm,1)      ;elw=dm_zeros(1,jm,1)      ;
    tbe=dm_zeros(jm,kb,1)     ;tbw=dm_zeros(jm,kb,1)     ;
    sbe=dm_zeros(jm,kb,1)     ;sbw=dm_zeros(jm,kb,1)     ;
    tbn=dm_zeros(im,kb,1)     ;tbs=dm_zeros(im,kb,1)     ;
    sbn=dm_zeros(im,kb,1)     ;sbs=dm_zeros(im,kb,1)     ;
    rfe            =0.0       ;rfn            =0.0  ;
    rfs            =0.0;                         

    dx_3d=dm_zeros(im,jm,kb)  ;dy_3d=dm_zeros(im,jm,kb) ;
    cor_3d=dm_zeros(im,jm,kb) ;
    h_3d=dm_zeros(im,jm,kb)   ;art_3d=dm_zeros(im,jm,kb);
    aru_3d=dm_zeros(im,jm,kb) ;arv_3d=dm_zeros(im,jm,kb) ;
    fsm_3d=dm_zeros(im,jm,kb) ;dum_3d=dm_zeros(im,jm,kb);
    dvm_3d=dm_zeros(im,jm,kb) ;dt_3d=dm_zeros(im,jm,kb);

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2D matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    aam2d=dm_zeros(im,jm,1)   ; advua=dm_zeros(im,jm,1)   ; 
    ady2d=dm_zeros(im,jm,1)   ; art=dm_zeros(im,jm,1)     ; 
    cbc=dm_zeros(im,jm,1)     ; cor=dm_zeros(im,jm,1)     ; 
    dry2d=dm_zeros(im,jm,1)   ; dt=dm_zeros(im,jm,1)      ; 
    dx=dm_zeros(im,jm,1)      ; dy=dm_zeros(im,jm,1)      ; 
    east_u=dm_zeros(im,jm,1)  ; east_v=dm_zeros(im,jm,1)  ; 
    egf=dm_zeros(im,jm,1)     ; el=dm_zeros(im,jm,1)      ; 
    et=dm_zeros(im,jm,1)      ; etb=dm_zeros(im,jm,1)     ; 
    fluxva=dm_zeros(im,jm,1)  ; fsm=dm_zeros(im,jm,1)     ; 
    north_e=dm_zeros(im,jm,1) ; north_u=dm_zeros(im,jm,1) ; 
    rot=dm_zeros(im,jm,1)     ; ssurf=dm_zeros(im,jm,1)   ; 
    tps=dm_zeros(im,jm,1)     ; tsurf=dm_zeros(im,jm,1)   ; 
    uab=dm_zeros(im,jm,1)     ; uaf=dm_zeros(im,jm,1)     ; 
    va=dm_zeros(im,jm,1)      ; vab=dm_zeros(im,jm,1)     ; 
    vtb=dm_zeros(im,jm,1)     ; vtf=dm_zeros(im,jm,1)     ; 
    wubot=dm_zeros(im,jm,1)   ; wusurf=dm_zeros(im,jm,1)  ; 
    wubot1=dm_zeros(im,jm,1)  ; wvbot1=dm_zeros(im,jm,1)  ; 

    advva=dm_zeros(im,jm,1)   ; adx2d=dm_zeros(im,jm,1)   ; 
    aru=dm_zeros(im,jm,1)     ; arv=dm_zeros(im,jm,1)     ; 
    d=dm_zeros(im,jm,1)       ; drx2d=dm_zeros(im,jm,1)   ; 
    dum=dm_zeros(im,jm,1)     ; dvm=dm_zeros(im,jm,1)     ; 
    east_c=dm_zeros(im,jm,1)  ; east_e=dm_zeros(im,jm,1)  ; 
    e_atmos=dm_zeros(im,jm,1) ; egb=dm_zeros(im,jm,1)     ; 
    elb=dm_zeros(im,jm,1)     ; elf=dm_zeros(im,jm,1)     ; 
    etf=dm_zeros(im,jm,1)     ; fluxua=dm_zeros(im,jm,1)  ; 
    h=dm_zeros(im,jm,1)       ; north_c=dm_zeros(im,jm,1) ; 
    north_v=dm_zeros(im,jm,1) ; psi=dm_zeros(im,jm,1)     ; 
    swrad=dm_zeros(im,jm,1)   ; vfluxb=dm_zeros(im,jm,1)  ; 
    ua=dm_zeros(im,jm,1)      ; vfluxf=dm_zeros(im,jm,1)  ; 
    utb=dm_zeros(im,jm,1)     ; utf=dm_zeros(im,jm,1)     ; 
    vaf=dm_zeros(im,jm,1)     ;                             
    wssurf=dm_zeros(im,jm,1)  ; wtsurf=dm_zeros(im,jm,1)  ; 
    wvbot=dm_zeros(im,jm,1)   ; wvsurf=dm_zeros(im,jm,1)  ;   


    !*************************
    !3D matrices
    !*************************
    aam=dm_zeros(im,jm,kb)  ; advx=dm_zeros(im,jm,kb) ;
    advy=dm_zeros(im,jm,kb) ; a=dm_zeros(im,jm,kb)    ; 
    c=dm_zeros(im,jm,kb)    ; drhox=dm_zeros(im,jm,kb);
    drhoy=dm_zeros(im,jm,kb); dtef=dm_zeros(im,jm,kb) ; 
    ee=dm_zeros(im,jm,kb)   ; gg=dm_zeros(im,jm,kb)   ;
    kh=dm_zeros(im,jm,kb)   ; km=dm_zeros(im,jm,kb)   ; 
    kq=dm_zeros(im,jm,kb)   ; l=dm_zeros(im,jm,kb)    ;
    q2b=dm_zeros(im,jm,kb)  ; q2=dm_zeros(im,jm,kb)   ; 
    q2lb=dm_zeros(im,jm,kb) ; q2l=dm_zeros(im,jm,kb)  ;
    rho=dm_zeros(im,jm,kb)  ; rmean=dm_zeros(im,jm,kb);
    sb=dm_zeros(im,jm,kb)   ; sclim=dm_zeros(im,jm,kb);
    s=dm_zeros(im,jm,kb)    ; tb=dm_zeros(im,jm,kb)   ; 
    tclim=dm_zeros(im,jm,kb); t=dm_zeros(im,jm,kb)    ;
    ub=dm_zeros(im,jm,kb)   ; uf=dm_zeros(im,jm,kb)   ; 
    u=dm_zeros(im,jm,kb)    ; vb=dm_zeros(im,jm,kb)   ;
    vf=dm_zeros(im,jm,kb)   ; v=dm_zeros(im,jm,kb)    ; 
    w=dm_zeros(im,jm,kb)    ; zflux=dm_zeros(im,jm,kb);


    ele=dm_zeros(1,jm,1)      ; eln=dm_zeros(1,im,1)      ;
    els=dm_zeros(1,im,1)      ; elw=dm_zeros(1,jm,1)      ;
    sbe=dm_zeros(jm,kb,1)     ; sbn=dm_zeros(im,kb,1)     ;
    sbs=dm_zeros(im,kb,1)     ; sbw=dm_zeros(jm,kb,1)     ; 
    tbe=dm_zeros(jm,kb,1)     ; tbn=dm_zeros(im,kb,1)     ;
    tbs=dm_zeros(im,kb,1)     ; tbw=dm_zeros(jm,kb,1)     ; 
    uabe=dm_zeros(1,jm,1)     ; uabw=dm_zeros(1,jm,1)     ;
    ube=dm_zeros(jm,kb,1)     ; ubw=dm_zeros(jm,kb,1)     ; 
    vabn=dm_zeros(1,im,1)     ; vabs=dm_zeros(1,im,1)     ;
    vbn=dm_zeros(im,kb,1)     ; vbs=dm_zeros(im,kb,1)     ; 

    d_3d=dm_zeros(im,jm,kb)   ; 

  end subroutine InitGrid

  subroutine FinalizeGrid()
    integer :: ierr
    
    call dm_destroy(dt, ierr);     call dm_destroy(tbe, ierr);
    call dm_destroy(dx_3d, ierr);
    call dm_destroy(dy_3d, ierr); call dm_destroy(cor_3d, ierr);
    call dm_destroy(h_3d, ierr); call dm_destroy(art_3d, ierr);
    call dm_destroy(aru_3d, ierr); call dm_destroy(arv_3d, ierr);
    
    call dm_destroy(fsm_3d, ierr); call dm_destroy(dum_3d, ierr);
    call dm_destroy(dvm_3d, ierr); call dm_destroy(dt_3d, ierr);
    call dm_destroy(d_3d, ierr);  

    call dm_destroy(aam2d, ierr); call dm_destroy(advua,ierr);
    call dm_destroy(advva, ierr); call dm_destroy(adx2d,ierr);
    call dm_destroy(ady2d, ierr); call dm_destroy(art, ierr);
    call dm_destroy(aru, ierr); call dm_destroy(arv, ierr);
    call dm_destroy(cbc, ierr); call dm_destroy(cor, ierr);
    call dm_destroy(d, ierr); call dm_destroy(drx2d, ierr);
    
    call dm_destroy(dry2d,ierr); call dm_destroy(dum,ierr);
    call dm_destroy(dvm,ierr); call dm_destroy(dx, ierr);
    call dm_destroy(dy, ierr); call dm_destroy(east_c,ierr);
    call dm_destroy(east_e,ierr); 
    call dm_destroy(east_u, ierr); call dm_destroy(east_v, ierr);
    call dm_destroy(e_atmos, ierr); call dm_destroy(egb, ierr);
    
    call dm_destroy(egf, ierr); call dm_destroy(el, ierr);
    call dm_destroy(elb, ierr); call dm_destroy(elf, ierr);
    call dm_destroy(et, ierr);  call dm_destroy(etb, ierr);
    call dm_destroy(etf, ierr);  call dm_destroy(fluxua, ierr);
    call dm_destroy(fluxva, ierr);  call dm_destroy(fsm, ierr);
    call dm_destroy(h,ierr);  call dm_destroy(north_c, ierr);
    
    call dm_destroy(north_e, ierr); call dm_destroy(north_u, ierr);
    call dm_destroy(north_v, ierr); call dm_destroy(psi, ierr);
    
    call dm_destroy(rot, ierr); call dm_destroy(ssurf, ierr);
    call dm_destroy(swrad, ierr);  call dm_destroy(vfluxb, ierr);
    call dm_destroy(tps, ierr);     call dm_destroy(tsurf, ierr);
    call dm_destroy(ua, ierr);   call dm_destroy(vfluxf, ierr);
    call dm_destroy(uab, ierr);  call dm_destroy(uaf, ierr);
    call dm_destroy(utb, ierr);  call dm_destroy(utf, ierr);     
    call dm_destroy(va, ierr);      call dm_destroy(vab, ierr);
    call dm_destroy(vaf, ierr);    
    call dm_destroy(vtb, ierr);  call dm_destroy(vtf, ierr);
    call dm_destroy(wssurf, ierr);  call dm_destroy(wtsurf, ierr); 
    call dm_destroy(wubot, ierr);   call dm_destroy(wusurf, ierr);
    call dm_destroy(wvbot, ierr);  call dm_destroy(wvsurf, ierr); 
    call dm_destroy(wubot1, ierr);  call dm_destroy(wvbot1, ierr); 

    call dm_destroy(aam, ierr); call dm_destroy(advx, ierr);
    call dm_destroy(advy, ierr); call dm_destroy(a, ierr);
    call dm_destroy(c, ierr);   call dm_destroy(drhox, ierr);
    call dm_destroy(drhoy, ierr); call dm_destroy(dtef, ierr); 
    call dm_destroy(ee, ierr);  call dm_destroy(gg, ierr);
    call dm_destroy(kh, ierr);  call dm_destroy(km, ierr);    
    call dm_destroy(kq, ierr);   call dm_destroy(l, ierr);
    call dm_destroy(q2b, ierr);  call dm_destroy(q2, ierr);    
    call dm_destroy(q2lb, ierr); call dm_destroy(q2l, ierr);
    call dm_destroy(rho, ierr);  call dm_destroy(rmean,ierr);
    call dm_destroy(sb, ierr); call dm_destroy(sclim,ierr);
    call dm_destroy(s, ierr);  call dm_destroy(tb, ierr);    
    call dm_destroy(tclim, ierr); call dm_destroy(t, ierr);
    call dm_destroy(ub, ierr);  call dm_destroy(uf, ierr);    
    call dm_destroy(u, ierr); call dm_destroy(vb, ierr);
    call dm_destroy(vf, ierr); call dm_destroy(v, ierr);    
    call dm_destroy(w, ierr); call dm_destroy(zflux, ierr);

    call dm_destroy(ele, ierr); call dm_destroy(eln, ierr);
    call dm_destroy(els, ierr); call dm_destroy(elw, ierr);       
    call dm_destroy(sbe,ierr);  call dm_destroy(sbn,ierr);
    call dm_destroy(sbs, ierr); call dm_destroy(sbw,ierr);      
    call dm_destroy(tbn, ierr); call dm_destroy(tbs, ierr);
    call dm_destroy(tbw, ierr);       
    call dm_destroy(uabe, ierr); call dm_destroy(uabw,ierr);
    
    call dm_destroy(ube, ierr);   call dm_destroy(ubw, ierr);    
    call dm_destroy(vabn, ierr);  call dm_destroy(vabs, ierr);
    call dm_destroy(vbn, ierr);   call dm_destroy(vbs, ierr);    

    call dm_destroy(z, ierr); call dm_destroy(zz, ierr);
    call dm_destroy(dz, ierr); call dm_destroy(dzz,ierr);
    call dm_destroy(z_3d,ierr); call dm_destroy(zz_3d,ierr);
    call dm_destroy(dz_3d,ierr); call dm_destroy(dzz_3d,ierr);

    call dm_destroy(curv, ierr)
    call dm_destroy(xflux, ierr)
    call dm_destroy(yflux, ierr)
    call dm_destroy(curv2d, ierr)
    call dm_destroy(etb_3d, ierr)
    call dm_destroy(egf_3d, ierr)
    call dm_destroy(egb_3d, ierr)
    call dm_destroy(e_atmos_3d, ierr)
    call dm_destroy(etf_3d, ierr)

    
  end subroutine FinalizeGrid
end module grid
