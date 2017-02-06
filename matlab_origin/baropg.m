function [rho,drhox,drhoy] = baropg(rho,drhox,drhoy,...
                             im,jm,imm1,jmm1,kb,kbm1,grav,...
                             zz,dt,dum,dvm,ramp,rmean,dx,dy)

% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
% *                                                                    *
% **********************************************************************
%     
      rho = rho - rmean;
%
%     Calculate x-component of baroclinic pressure gradient:
%
      for j=2:jmm1
        for i=2:imm1
          drhox(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))  ...
                  *(rho(i,j,1)-rho(i-1,j,1));
        end
      end
%
      for k=2:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhox(i,j,k)=drhox(i,j,k-1)     ...
                    +grav*.25e0*(zz(k-1)-zz(k))     ...
                      *(dt(i,j)+dt(i-1,j))     ...
                      *(rho(i,j,k)-rho(i-1,j,k)     ... 
                       +rho(i,j,k-1)-rho(i-1,j,k-1))     ...   
                 +grav*.25e0*(zz(k-1)+zz(k))     ...
                      *(dt(i,j)-dt(i-1,j))     ... 
                     *(rho(i,j,k)+rho(i-1,j,k)     ...
                        -rho(i,j,k-1)-rho(i-1,j,k-1));
          end
        end
      end
%
      for k=1:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))     ... 
                       *drhox(i,j,k)*dum(i,j)     ...   
                     *(dy(i,j)+dy(i-1,j));
          end
        end
      end
%
%     Calculate y-component of baroclinic pressure gradient:
%
      for j=2:jmm1
        for i=2:imm1
          drhoy(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))     ...  
                *(rho(i,j,1)-rho(i,j-1,1));
        end
      end
%
      for k=2:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)     ... 
                   +grav*.25e0*(zz(k-1)-zz(k))     ...
                      *(dt(i,j)+dt(i,j-1))     ...  
                    *(rho(i,j,k)-rho(i,j-1,k)     ... 
                       +rho(i,j,k-1)-rho(i,j-1,k-1))     ...
                    +grav*.25e0*(zz(k-1)+zz(k))     ...
                      *(dt(i,j)-dt(i,j-1))     ...
                      *(rho(i,j,k)+rho(i,j-1,k)     ...
                        -rho(i,j,k-1)-rho(i,j-1,k-1));
          end
        end
      end
%
      for k=1:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))     ...
                        *drhoy(i,j,k)*dvm(i,j)     ...
                        *(dx(i,j)+dx(i,j-1));
          end
        end
      end
%
      for k=1:kb
        for j=2:jmm1
          for i=2:imm1
            drhox(i,j,k)=ramp*drhox(i,j,k);
            drhoy(i,j,k)=ramp*drhoy(i,j,k);
          end
        end
      end
%
      for k=1:kb
        for j=1:jm
          for i=1:im
            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k);
          end
        end
      end
return
end