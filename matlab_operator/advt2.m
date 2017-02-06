function [fb,f,fclim,ff,xflux,yflux,nitera,sw,...
    zflux] = advt2(fb,f,fclim,ff,xflux,yflux,nitera,sw,...
                                                 zflux,...
                                                 im,jm,kb,imm1,jmm1,kbm1,dti2,...
                                                 etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz)



% **********************************************************************
% *                                                                    *
% * FUN%TION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is a first-order upstream scheme, which        *
% *                reduces implicit diffusion using the Smolarkiewicz  *
% *                iterative upstream scheme with an antidiffusive     *
% *                velocity.                                           *
% *                                                                    *
% *                It is based on the subroutines of Gianmaria Sannino *
% *                (Inter-university %omputing %onsortium, Rome, Italy)*
% *                and Vincenzo Artale (Italian National Agency for    *
% *                New Technology and Environment, Rome, Italy),       *
% *                forwnloaded from the POM FTP site on 1 Nov. 2001.    *
% *                The calculations have been simplified by removing   *
% *                the shock switch option. It should be noted that    *
% *                this implementation fores not include cross-terms    *
% *                which are in the original formulation.              *
% *                                                                    *
% *                fb,f,fclim,ff . as used in subroutine advt1         *
% *                xflux,yflux ... working arrays used to save memory  *
% *                nitera ........ number of iterations. This should   *
% *                                be in the range 1 - 4. 1 is         *
% *                                standard upstream differencing;     *
% *                                3 adds 50% %PU time to POM.         *
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
% *                  %omputational Physics, 54, 325-362, 1984.         *
% *                                                                    *
% **********************************************************************
%
%
%
%
%     calculate horizontal mass fluxes:
%
xmassflux=zeros(im,jm,kb);
ymassflux=zeros(im,jm,kb);


%
for k=1:kbm1
    for j=2:jmm1
        for i=2:im
            xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j)) ...
                *(dt(i-1,j)+dt(i,j))*u(i,j,k);
        end
    end
    %
    for j=2:jm
        for i=2:imm1
            ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j))     ...
                *(dt(i,j-1)+dt(i,j))*v(i,j,k);
        end
    end
end
%
for j=1:jm
    for i=1:im
        fb(i,j,kb)=fb(i,j,kbm1);
    end
end
%
eta=etb;
zwflux=w;
fbmem=fb;
%
%     Start Smolarkiewicz scheme:
%
for itera=1:nitera
    %
    %     Upwind advection scheme:
    %
    for k=1:kbm1
        for j=2:jm
            for i=2:im
                xflux(i,j,k)=0.5e0     ...
                    *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))     ...
                    *fbmem(i-1,j,k)+     ...
                    (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))     ...
                    *fbmem(i,j,k));
                %
                yflux(i,j,k)=0.5e0     ...
                    *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))     ...
                    *fbmem(i,j-1,k)+     ...
                    (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))     ...
                    *fbmem(i,j,k));
            end
        end
    end
    %
    for j=2:jmm1
        for i=2:imm1
            zflux(i,j,1)=0.e0;
            if(itera==1)
                zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j);
            end
            zflux(i,j,kb)=0.e0;
        end
    end
    %
    for k=2:kbm1
        for j=2:jmm1
            for i=2:imm1
                zflux(i,j,k)=0.5e0     ...
                    *((zwflux(i,j,k)+abs(zwflux(i,j,k)))     ...
                    *fbmem(i,j,k)+     ...
                    (zwflux(i,j,k)-abs(zwflux(i,j,k)))     ...
                    *fbmem(i,j,k-1));
                zflux(i,j,k)=zflux(i,j,k)*art(i,j);
            end
        end
    end
    %
    %     Add net advective fluxes and step forward in time:
    %
    for k=1:kbm1
        for j=2:jmm1
            for i=2:imm1
                ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)     ...
                    +yflux(i,j+1,k)-yflux(i,j,k)     ...
                    +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k);
                ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)     ...
                    -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j));
            end
        end
    end
    %
    %     %alculate antidiffusion velocity:
    %
  
   [xmassflux,ymassflux,zwflux,ff,sw] = smol_adif(xmassflux,ymassflux,zwflux,ff,sw,...
                                                        im,jm,kb,imm1,jmm1,kbm1,dti2,fsm,aru,arv,dzz,dt);
    %
    eta=etf;
    fbmem=ff;
    %
    %     End of Smolarkiewicz scheme
    %
end
%
%     Add horizontal diffusive fluxes:
%

fb=fb-fclim;
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k));
            ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k));
        end
    end
end
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni     ...
                *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)     ...
                *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j));
            yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni     ...
                *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)     ...
                *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1));
        end
    end
end
%

fb=fb+fclim;
%
%     Add net horizontal fluxes and step forward in time:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)     ...
                +yflux(i,j+1,k)-yflux(i,j,k))     ...
                /((h(i,j)+etf(i,j))*art(i,j));
        end
    end
end

end

%
% **********************************************************************
% *                                                                    *
