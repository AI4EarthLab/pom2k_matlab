function [drhox,drhoy] = baropg_mcc(rho,rmean,d,dt,ramp)
% calculate  baroclinic pressure gradient 4th order correction terms, following McCalpin
global im jm kb imm1 jmm1 kbm1 dum dvm dum_3d dvm_3d dx_3d dx dy dy_3d dz_3d grav zz_3d gs;
dt_3d=repmat(dt,1,1,kb);
% drhoy=zeros(im,jm,kb);      drhox=zeros(im,jm,kb);
% drho=create_field(zeros(im,jm,kb),gs,3);
    rho=rho-rmean;


% convert a 2nd order matrices to special 4th order 
% special 4th order case
%       call order2d_mpi(d,d4th,im,jm)
% d4th(2:im+1,2:jm+1)=d;          d4th(1,2:jm+1)=d(im-2,:);    d4th(2:im+1,1)=d(:,jm-2);
%       call order3d_mpi(rho,rho4th,im,jm,kb)
% d4th=d; rho4th=rho;
% ! compute terms correct to 4th order
% ! compute DRHO, RHOU, DDX and D4
%       for j=1:jm
%         for i=2:im
%           for k=1:kbm1
%             drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j);
%             rhou(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j);
%           end
%           ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j);
%           d4(i,j)=.5*(d(i,j)+d(i-1,j))*dum(i,j);
%         end
%       end
drho=   DXB(rho).*AXB(dx_3d).*dum_3d;
rhou=   AXB(rho).*dum_3d;
ddx=    DXB(d).*AXB(dx).*dum;
d4=     AXB(d).*dum;
drho=create_field(drho.data,gs,3);     
      
%       if(n_west == -1)
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
        
%! calculate x-component of baroclinic pressure gradient
%       for j=2:jmm1
%         for i=2:imm1
%           drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1);
%         end
%       end
% 
%       for k=2:kbm1
%         for j=2:jmm1
%           for i=2:imm1
%             drhox(i,j,k)=drhox(i,j,k-1)+grav*0.5e0*dzz(k-1)*d4(i,j)*(drho(i,j,k-1)+drho(i,j,k)) ...
%                         +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)*(rhou(i,j,k)-rhou(i,j,k-1));
%           end
%         end
%       end

% drhox=ramp * grav * AXB(dt_3d) .* SUMZ(-DZB(zz_3d).* AXB(dt_3d) .* DXB(AZB(rho-rmean)) + AZB(zz_3d) .* DXB(dt_3d).* DZB(AXB(rho-rmean)).*AZB(dz_3d)).* dum_3d;
% drho.pos=3;
drhox=AXB(dt_3d) .* ramp .* grav .* SUMZ(-DZB(zz_3d).* repmat(d4,1,1,kb) .* AZB(drho) +  AZB(zz_3d) .* repmat(ddx,1,1,kb) .* DZB(rhou).*AZB(dz_3d)).* dum_3d .* AXB(dy_3d);

%       for k=1:kbm1
%         for j=2:jmm1
%           for i=2:imm1
%             drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))*drhox(i,j,k)*dum(i,j)*(dy(i,j)+dy(i-1,j));
%           end 
%         end
%       end

% ! compute terms correct to 4th order
% ddx=zeros(im,jm);   d4=zeros(im,jm);    rhou=zeros(im,jm,kb);   drho=zeros(im,jm,kb);
      
% ! compute DRHO, RHOU, DDX and D4
%       for j=2:jm
%         for i=1:im
%           for k=1:kbm1
%             drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j);
%             rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j);
%           end
%           ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j);
%           d4(i,j)=.5*(d(i,j)+d(i,j-1))*dvm(i,j);
%         end
%       end
    drho=   DYB(rho).*AYB(dy_3d).*dvm_3d;
    rhou=   AYB(rho).*dvm_3d;
    ddx=    DYB(d).*AYB(dy).*dvm;
    d4=     AYB(d).*dvm;
      
    drho=create_field(drho.data,gs,3);
%       if(n_south == -1)
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

% ! calculate y-component of baroclinic pressure gradient
%       for j=2:jmm1
%         for i=2:imm1
%           drhoy(i,j,1)=grav*(-zz(1)).*d4(i,j).*drho(i,j,1);
%         end
%       end
% 
%       for k=2:kbm1
%         for j=2:jmm1
%           for i=2:imm1
%             drhoy(i,j,k)=drhoy(i,j,k-1) +grav*0.5e0*dzz(k-1)*d4(i,j)*(drho(i,j,k-1)+drho(i,j,k))+grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)*(rhou(i,j,k)-rhou(i,j,k-1));
%           end
%         end
%       end

% drho.pos=3;
drhoy=AYB(dt_3d) .* ramp .* grav .* SUMZ(-DZB(zz_3d).* repmat(d4,1,1,kb) .* AZB(drho) +  AZB(zz_3d) .* repmat(ddx,1,1,kb) .* DZB(rhou).*AZB(dz_3d)).* dvm_3d .* AYB(dx_3d);
            
%       for k=1:kbm1
%         for j=2:jmm1
%           for i=2:imm1
%             drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1)) *drhoy(i,j,k)*dvm(i,j) *(dx(i,j)+dx(i,j-1));
%           end
%         end
%       end
    
      end