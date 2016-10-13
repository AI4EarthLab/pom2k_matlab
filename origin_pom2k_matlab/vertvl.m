function [xflux,yflux,...
    w]=vertvl (xflux,yflux,...
    w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  %alculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
%
%
%     Reestablish boundary conditions:
%
      for k=1:kbm1
        for j=2:jm
          for i=2:im
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)) ...
                    *(dt(i,j)+dt(i-1,j))*u(i,j,k);
          end
        end
      end
%
      for k=1:kbm1
        for j=2:jm
          for i=2:im
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1))     ... 
                   *(dt(i,j)+dt(i,j-1))*v(i,j,k);
          end
        end
      end
%
%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
%
        for j=2:jmm1
          for i=2:imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j));
          end
        end
%
      for k=1:kbm1
        for j=2:jmm1
          for i=2:imm1
            w(i,j,k+1)=w(i,j,k)     ...
                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)     ...  
                      +yflux(i,j+1,k)-yflux(i,j,k))     ... 
                       /(dx(i,j)*dy(i,j))     ... 
                       +(etf(i,j)-etb(i,j))/dti2);
          end
        end
      end
%

return
end
%
%   
%
%-----------------------------------------------------------------------

