function [xflux,yflux,...
    w]=vertvl (xflux,yflux,...
    w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  %alculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************

xflux = zeros(im,jm,kb);
yflux = zeros(im,jm,kb);
      for k=1:kbm1
          xflux(:,:,k) = AXB1_XY(dy) .* AXB1_XY(dt) .* u(:,:,k);
      end
%
      for k=1:kbm1
            yflux(:,:,k)=AYB1_XY(dx) .* AYB1_XY(dt) .* v(:,:,k);
      end
%
%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
%
w = zeros(im,jm,kb);
w(:,:,1)=0.5*(vfluxb+vfluxf);

%etf-etb在边界处不能有值
      for k=1:kbm1
            w(:,:,k+1)=w(:,:,k) + dz(k) ...
                .* ( ( DXF2_XY(xflux(:,:,k))+DYF2_XY(yflux(:,:,k)) )...
                ./(dx.*dy)+OP_L_XY*(etf-etb)/dti2 );
      end
%

return
end



