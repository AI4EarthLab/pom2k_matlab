function[xmassflux,ymassflux,zwflux,ff] = smol_adif(xmassflux,ymassflux,zwflux,ff,sw,aru,arv,dt)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  %calculates the antidiffusive velocity used to       *
% *                reduce the numerical diffusion associated with the  *
% *                upstream differencing scheme.                       *
% *                                                                    *
% *                This is based on a subroutine of Gianmaria Sannino  *
% *                (Inter-university computing consortium, Rome, Italy)*
% *                and Vincenzo Artale (Italian National Agency for    *
% *                New Technology and Environment, Rome, Italy),       *
% *                forwnloaded from the POM FTP site on 1 Nov. 2001.    *
% *                The calculations have been simplified by removing   *
% *                the shock switch option.                            *
% *                                                                    *
% **********************************************************************
global im jm kb imm1 jmm1 kbm1 dti2 fsm_3d dzz;
	value_min=1.0e-9;   epsilon=1.0e-14;
	ff=ff.*fsm_3d;      dt_3d=repmat(dt,1,1,kb);
%
%     Recalculate mass fluxes with antidiffusion velocity:
%
      for k=1:kbm1
        for j=2:jmm1
          for i=2:im
            if(ff.data(i,j,k)<value_min || ff.data(i-1,j,k)<value_min)
              xmassflux(i,j,k)=0.e0;
            else
              udx=abs(xmassflux(i,j,k));
              u2dt=dti2*xmassflux(i,j,k).^2 * 2.e0 ./(aru(i,j)*(dt(i-1,j)+dt(i,j)));
              mol=(ff(i,j,k)-ff(i-1,j,k))/(ff(i-1,j,k)+ff(i,j,k)+epsilon);
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw;
              if(abs(udx.data)<abs(u2dt.data)) 
                  xmassflux(i,j,k)=0.e0;
              end 
            end
          end
        end
      end

      for k=1:kbm1
        for j=2:jm
          for i=2:imm1
            if(ff.data(i,j,k)<value_min|| ff.data(i,j-1,k)<value_min) 
              ymassflux(i,j,k)=0.e0;
            else
             vdy=abs(ymassflux(i,j,k));
             v2dt=dti2*ymassflux(i,j,k).^2 *2.e0 /(arv(i,j)*(dt(i,j-1)+dt(i,j)));
             mol=(ff(i,j,k)-ff(i,j-1,k)) /(ff(i,j-1,k)+ff(i,j,k)+epsilon);
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw;
             if(abs(vdy.data)<abs(v2dt.data))
				 ymassflux(i,j,k)=0.e0;
             end
            end
          end
        end
      end
%
      for k=2:kbm1
        for j=2:jmm1
          for i=2:imm1
            if(ff.data(i,j,k)<value_min|| ff.data(i,j,k-1)<value_min)
              zwflux(i,j,k)=0.e0;
            else
              wdz=abs(zwflux(i,j,k));
              w2dt=dti2*zwflux(i,j,k).^2/(dzz(k-1)*dt(i,j));
              mol=(ff(i,j,k-1)-ff(i,j,k)) /(ff(i,j,k)+ff(i,j,k-1)+epsilon);
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw;
              if(abs(wdz.data)<abs(w2dt.data))
				zwflux(i,j,k)=0.e0;
              end 
            end 
          end
        end
      end
      
end
