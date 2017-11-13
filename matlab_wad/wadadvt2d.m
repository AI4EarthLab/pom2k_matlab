      function [ff]=wadadvt2d(fb,f,ff,nitera,im,jm,kb,sw,dx,dy,u,v,art,aru,arv, ...
         fsm,dum,dvm,aam,dti2) 
%{
%     2-d version of Smolarkiewicz from /home/lyo/pom/pom2k/           !
%     pom2k.f's subroutine advt2:                                      !                                                               !
%     This version gets rid of common block, i.e. the routine is now   !
%     self-contained, except that im & jm need to be specified below   !
%     if "include 'grid'" is NOT used; also, 2-d calculation, i.e,     !
%     solve for  D:                                                    !
%                                                                      !
%     d(D)/dt + d(UD)/dx + d(VD)/dy = Diffusion_of_(D)                 !
%                                                                      !
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is a first-order upstream scheme, which        *
% *                reduces implicit diffusion using the Smolarkiewicz  *
% *                iterative upstream scheme with an antidiffusive     *
% *                velocity.                                           *
% *                                                                    *
% *                It is based on the subroutines of Gianmaria Sannino *
% *                (Inter-university Computing Consortium, Rome, Italy)*
% *                and Vincenzo Artale (Italian National Agency for    *
% *                New Technology and Environment, Rome, Italy),       *
% *                downloaded from the POM FTP site on 1 Nov. 2001.    *
% *                The calculations have been simplified by removing   *
% *                the shock switch option. It should be noted that    *
% *                this implementation does not include cross-terms    *
% *                which are in the original formulation.              *
% *                                                                    *
% *                fb,f,fclim,ff . as used in subroutine advt1         *
% *                xflux,yflux ... working arrays used to save memory  *
% *                nitera ........ number of iterations. This should   *
% *                                be in the range 1 - 4. 1 is         *
% *                                standard upstream differencing;     *
% *                                3 adds 50% CPU time to POM.         *
% *                sw ............ smoothing parameter. This should    *
% *                                preferably be 1, but 0 < sw < 1     *
% *                                gives smoother solutions with less  *
% *                                overshoot when nitera > 1.          *
% *                                                                    *
% *                Reference:                                          *
% *                                                                    *
% *                Smolarkiewicz, P.K.; A fully multidimensional       *
% *                  positive definite advection transport algorithm   *
% *                  with small implicit diffusion, Journal of         *
% *                  Computational Physics, 54, 325-362, 1984.         *
% *                                                                    *
% **********************************************************************
%}

      fbmem=zeros(im,jm); 
      xflux=fbmem     ; yflux=fbmem ;
      xmassflux=fbmem ; ymassflux=fbmem ;
      imm1=im-1; jmm1=jm-1;
      xmassflux(2:im  ,2:jmm1)=0.5e0*(dy(1:imm1,2:jmm1)+dy(2:im  ,2:jmm1))*u(2:im  ,2:jmm1);
      xmassflux(2:imm1,2:jm  )=0.5e0*(dx(2:imm1,1:jmm1)+dy(2:imm1,2:jm  ))*u(2:imm1,2:jm  );
      fbmem=fb;

%     Start Smolarkiewicz scheme:

      for itera=1:nitera

%     Upwind advection scheme:

          for j=2:jm
            for i=2:im
              vt=xmassflux(i,j)
              xflux(i,j)=0.5*((vt+abs (vt))*fbmem(i-1,j)+(vt-abs (vt))*fbmem(i,j);
              vt=ymassflux(i,j)
              yflux(i,j)=0.5*((vt+abs (vt))*fbmem(i,j-1)+(vt-abs (vt))*fbmem(i,j);
            end 
          end 
%     Add net advective fluxes and step forward in time:

          for j=2:jmm1
            for i=2:imm1
              ff(i,j)=xflux(i+1,j)-xflux(i,j)+yflux(i,j+1)-yflux(i,j);
              ff(i,j)=(fbmem(i,j)*art(i,j)-dti2*ff(i,j))/(art(i,j));
            end
          end 

%     Calculate antidiffusion velocity:

			[xmassflux,ymassflux,ff]=wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2,im,km)
      fbmem =fb;

%     End of Smolarkiewicz scheme

      end 

%     Add horizontal diffusive fluxes:

      for j=2:jm
        for i=2:im
          xmassflux(i,j)=0.5e0*(aam(i,j)+aam(i-1,j));
          ymassflux(i,j)=0.5e0*(aam(i,j)+aam(i,j-1));
        end
      end
        for j=2,jm
          for i=2,im
           xflux(i,j)=-xmassflux(i,j)           ...
                *(fb(i,j)-fb(i-1,j))*dum(i,j)   ...
                *(dy(i,j)+dy(i-1,j))/(dx(i,j)+dx(i-1,j));
           yflux(i,j)=-ymassflux(i,j)           ...
                *(fb(i,j)-fb(i,j-1))*dvm(i,j)   ...
                *(dx(i,j)+dx(i,j-1))/(dy(i,j)+dy(i,j-1));
          end
        end

%     Add net horizontal fluxes and step forward in time:

        for j=2:jmm1
          for i=2:imm1
            ff(i,j)=ff(i,j)-dti2*(xflux(i+1,j)-xflux(i,j) ...
          +yflux(i,j+1)-yflux(i,j))/(art(i,j));
          end
        end
end 
