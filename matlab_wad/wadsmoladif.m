 function [xmassflux,ymassflux,ff]=wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2,im,km)
%                                                                       %  
%  **********************************************************************
%  *                                                                    *
%  * FUNCTION    :  Calculates the antidiffusive velocity used to       *
%  *                reduce the numerical diffusion associated with the  *
%  *                upstream differencing scheme.                       *
%  *                                                                    *
%  *                This is based on a subroutine of Gianmaria Sannino  *
%  *                (Inter-university Computing Consortium, Rome, Italy)*
%  *                and Vincenzo Artale (Italian National Agency for    *
%  *                New Technology and Environment, Rome, Italy),       *
%  *                downloaded from the POM FTP site on 1 Nov. 2001.    *
%  *                The calculations have been simplified by removing   *
%  *                the shock switch option.                            *
%  *                                                                    *
%  **********************************************************************
%  
%      integer im,jm,kb  % lyo:% wad:kb not required but defined
                        %          in "grid" below
% 
%      PARAMETER (IM=65,JM=49)
%      PARAMETER (IM=131,JM=99)
%      include 'grid'
% 
%      real ff(im,jm)
%      real xmassflux(im,jm),ymassflux(im,jm)
%      real fsm(im,jm),aru(im,jm),arv(im,jm)
%      real sw,dti2
%      real mol,abs_1,abs_2
%      real value_min,epsilon
%      real udx,u2dt,vdy,v2dt,wdz,w2dt
%      integer i,j,k,imm1,jmm1
% 
      value_min=1.e-9;epsilon=1.0e-14;
% 
% 
      imm1=im-1; jmm1=jm-1;
% 
%      Apply temperature and salinity mask:
% 
        for i=1:im
          for j=1:jm
            ff(i,j)=ff(i,j)*fsm(i,j);
          end
        end
% 
%      Recalculate mass fluxes with antidiffusion velocity:
% 
        for j=2:jmm1
          for i=2:im
            if(ff(i,j)<value_min||ff(i-1,j)<value_min) 
              xmassflux(i,j)=0.e0;
            else
              udx=abs(xmassflux(i,j));
              u2dt=dti2*xmassflux(i,j)*xmassflux(i,j)/(aru(i,j));
              mol=(ff(i,j)-ff(i-1,j))/(ff(i-1,j)+ff(i,j)+epsilon);
              xmassflux(i,j)=(udx-u2dt)*mol*sw;
              abs_1=abs(udx);
              abs_2=abs(u2dt);
              if(abs_1<abs_2)
                xmassflux(i,j)=0.e0
              end
            end
          end
        end
% 
        for j=2:jm
          for i=2:imm1
            if(ff(i,j)<value_min.or.ff(i,j-1)<value_min) 
              ymassflux(i,j)=0.e0;
            else
             vdy=abs(ymassflux(i,j));
             v2dt=dti2*ymassflux(i,j)*ymassflux(i,j)/(arv(i,j));
             mol=(ff(i,j)-ff(i,j-1))/(ff(i,j-1)+ff(i,j)+epsilon);
             ymassflux(i,j)=(vdy-v2dt)*mol*sw;
             abs_1=abs(vdy);      abs_2=abs(v2dt);
             if(abs_1<abs_2)
             	 ymassflux(i,j)=0.e0;
             end
            end
          end
        end  
end 
