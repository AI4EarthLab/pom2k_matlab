
      function [a,c,ee,gg,tps,vf,wvbot] = profv(a,c,ee,gg,tps,vf,wvbot,...
                                          dvm,dz,dzz,im,jm,kb,imm1,jmm1,kbm1,kbm2,...
                                          km,cbc,ub,vb,umol,wvsurf,h,etf,dti2)
                                            
% **********************************************************************
%                                                                      *
% * FUNCTION    :  Solves for vertical diffusion of y-momentum using   *
% *                method described by Richmeyer and Morton.           *
% *                                                                    *
% *                See:                                                *
% *                                                                    *
% *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
% *                  Methods for Initial-Value Problems, 2nd edition,  *
% *                  Interscience, New York, 198-201.                  *
% *                                                                    *
% *                NOTE that wvsurf has the opposite sign to the wind  *
% *                speed.                                              *
% *                                                                    *
% **********************************************************************
%
%
	 dh=ones(im,jm);
%     The following section solves the equation:
%
%       dti2*(km*u')'-u=-ub
%
%
      for j=2:jm
        for i=2:im
          dh(i,j)=.5e0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1));
        end
      end
%
      for k=1:kb
        for j=2:jm
          for i=2:im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5e0;
          end
        end
      end
%
      for k=2:kbm1
        for j=1:jm
          for i=1:im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)     ...
                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j));
            c(i,j,k)=-dti2*(c(i,j,k)+umol)     ...
                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j));
          end
        end
      end
%
      for j=1:jm
        for i=1:im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0);
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))     ... 
              /(a(i,j,1)-1.e0);
        end
      end
%
      for k=2:kbm2
        for j=1:jm
          for i=1:im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0);
            ee(i,j,k)=a(i,j,k)*gg(i,j,k);
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k);
          end
        end
      end
%
      for j=2:jmm1
        for i=2:imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i,j-1))     ...
              *sqrt((.25e0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)     ...
                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))^2     ...
                    +vb(i,j,kbm1)^2);
     
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))     ... 
                 /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0     ...  
                  -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1));
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j);
        end
      end
%
      for k=2:kbm1
        ki=kb-k;
        for j=2:jmm1
          for i=2:imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j);
          end
        end
      end
%
      for j=2:jmm1
        for i=2:imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1);
        end
      end
%
      return
%
      end
%
