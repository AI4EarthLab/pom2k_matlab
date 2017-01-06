
program pom2k
  use dm_op 
  use dm
  use dm_type
  use constants
  
  implicit none

  integer :: ierr
  integer :: myrank, mysize
  
  type(Matrix) aam2d, advua, advva, adx2d, ady2d, art, aru, arv, cbc, cor, d, drx2d
  type(Matrix) dry2d   , dt      , dum     , dvm     
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
  type(Matrix) w    , zflux

  type(Matrix) ele      , eln      , els   , elw       
  type(Matrix) sbe     , sbn     , sbs     , sbw      
  type(Matrix) tbe     , tbn     , tbs     , tbw      
  type(Matrix) uabe     , uabw     , ube     , ubw    
  type(Matrix) vabn     , vabs     , vbn     , vbs    

  type(Matrix) d_3d, dt_3d
  
  call dm_init1(ierr)

  call dm_comm_rank(myrank,ierr)
  call dm_comm_size(mysize,ierr)

  call dm_init2(im, jm , kb,ierr)
  
  if(myrank==0) then 
     print *, "==============Input paramenters==========="
     print *, "im=",im,",jm=",jm,",kb=",kb
  endif

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



  call dm_finalize(ierr)
end program pom2k
