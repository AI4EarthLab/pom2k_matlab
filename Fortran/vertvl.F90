

subroutine vertvl(dti2)
  use dm_op
  use grid
  use input
  use dm
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

end subroutine 
