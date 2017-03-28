function [drhox,drhoy] = baropg_mcc(rho,rmean,d,ramp)
% calculate  baroclinic pressure gradient 4th order correction terms, following McCalpin
global im jm kb imm1 jmm1 kbm1

ddx=zeros(im,jm);   d4=zeros(im,jm);    rhou=zeros(im,jm,kb);   drho=zeros(im,jm,kb);
rho4th=zeros(im+1,jm+1,kb);             d4th=zeros(im+1,jm+1);  drhoy=zeros(im,jm,kb);
drhox=zeros(im,jm,kb);
%rho4th(0:im,0:jm,kb),d4th(0:im,0:jm)

    rho=rho-rmean;


% convert a 2nd order matrices to special 4th order 
% special 4th order case
%       call order2d_mpi(d,d4th,im,jm)
% d4th(2:im+1,2:jm+1)=d;          d4th(1,2:jm+1)=d(im-2,:);    d4th(2:im+1,1)=d(:,jm-2);
%       call order3d_mpi(rho,rho4th,im,jm,kb)

% ! compute terms correct to 4th order
% ! compute DRHO, RHOU, DDX and D4
      for j=1:jm
        for i=2:im
          for k=1:kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j);
            rhou(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j);
          end
          ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j);
          d4(i,j)=.5*(d(i,j)+d(i-1,j))*dum(i,j);
        end
      end

      if(n_west == -1)
        for j=1:jm
          for i=3:imm1
            for k=1:kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*(dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-2*(rho(i,j,k)-rho(i-1,j,k))+  ...
                         dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)));
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)* (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)));
            end
            ddx(i,j)=ddx(i,j)-(1./24.)*(dum(i+1,j)*(d(i+1,j)-d(i,j))-2*(d(i,j)-d(i-1,j))+ dum(i-1,j)*(d(i-1,j)-d(i-2,j)));
            d4(i,j)=d4(i,j)+(1./16.)* (dum(i+1,j)*(d(i,j)-d(i+1,j))+ dum(i-1,j)*(d(i-1,j)-d(i-2,j)));
          end
        end
      else
        for j=1:jm
          for i=2:imm1
            for k=1:kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*(dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))- 2*(rho(i,j,k)-rho(i-1,j,k))+ ...      ...
                         dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)));
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*(dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+ dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)));
            end
            ddx(i,j)=ddx(i,j)-(1./24.)* (dum(i+1,j)*(d(i+1,j)-d(i,j))- 2*(d(i,j)-d(i-1,j))+ dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)));
            d4(i,j)=d4(i,j)+(1./16.)*(dum(i+1,j)*(d(i,j)-d(i+1,j))+ dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)));
          end 
        end
      end

%! calculate x-component of baroclinic pressure gradient
      for j=2:jmm1
        for i=2:imm1
          drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1);
        end
      end

      for k=2:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhox(i,j,k)=drhox(i,j,k-1)+grav*0.5e0*dzz(k-1)*d4(i,j)*(drho(i,j,k-1)+drho(i,j,k)) ...
                        +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)*(rhou(i,j,k)-rhou(i,j,k-1));
          end
        end
      end

      for k=1:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))*drhox(i,j,k)*dum(i,j)*(dy(i,j)+dy(i-1,j));
          end 
        end
      end

% ! compute terms correct to 4th order
ddx=zeros(im,jm);   d4=zeros(im,jm);    rhou=zeros(im,jm,kb);   drho=zeros(im,jm,kb);
      
% ! compute DRHO, RHOU, DDX and D4
      for j=2:jm
        for i=1:im
          for k=1:kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j);
            rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j);
          end
          ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j);
          d4(i,j)=.5*(d(i,j)+d(i,j-1))*dvm(i,j);
        end
      end

      if(n_south == -1)
        for j=3:jmm1
          for i=1:im
            for k=1:kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*(dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))- ...
                         2*(rho(i,j,k)-rho(i,j-1,k))+dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)));
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*(dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+ ...
                         dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)));
            end
            ddx(i,j)=ddx(i,j)-(1./24)*(dvm(i,j+1)*(d(i,j+1)-d(i,j))-2*(d(i,j)-d(i,j-1))+dvm(i,j-1)*(d(i,j-1)-d(i,j-2)));
            d4(i,j)=d4(i,j)+(1./16.)*(dvm(i,j+1)*(d(i,j)-d(i,j+1))+dvm(i,j-1)*(d(i,j-1)-d(i,j-2)));
          end
        end 
      else
        for j=2:jmm1
          for i=1:im
            for k=1:kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*(dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))-2*(rho(i,j,k)-rho(i,j-1,k))+ ...
                         dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)));
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*(dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)));
            end
            ddx(i,j)=ddx(i,j)-(1./24)*(dvm(i,j+1)*(d(i,j+1)-d(i,j))-2*(d(i,j)-d(i,j-1))+ ...
                    dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)));
            d4(i,j)=d4(i,j)+(1./16.)*(dvm(i,j+1)*(d(i,j)-d(i,j+1))+dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)));
          end
        end
      end

% ! calculate y-component of baroclinic pressure gradient
      for j=2:jmm1
        for i=2:imm1
          drhoy(i,j,1)=grav*(-zz(1)).*d4(i,j).*drho(i,j,1);
        end
      end

      for k=2:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhoy(i,j,k)=drhoy(i,j,k-1) +grav*0.5e0*dzz(k-1)*d4(i,j)* ...
                (drho(i,j,k-1)+drho(i,j,k))+grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)*(rhou(i,j,k)-rhou(i,j,k-1));
          end
        end
      end

      for k=1:kbm1
        for j=2:jmm1
          for i=2:imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1)) *drhoy(i,j,k)*dvm(i,j) *(dx(i,j)+dx(i,j-1));
          end
        end
      end

    drhox=ramp.*drhox;
    drhoy=ramp.*drhoy;
    rho=rho+rmean;

      return
      end