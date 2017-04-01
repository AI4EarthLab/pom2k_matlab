function[xmassflux,ymassflux,zwflux,ff] = smol_adif(xmassflux,ymassflux,zwflux,ff,sw,dt)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  calculates the antidiffusive velocity used to       *
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
global kb dti2 fsm_3d dz_3d dx_3d dy_3d aru_3d arv_3d dzz_3d;
	value_min=1.0e-9;   epsilon=1.0e-14;
	ff=ff.*fsm_3d;      dt_3d=repmat(dt,1,1,kb);

%     Recalculate mass fluxes with antidiffusion velocity:
      xmassflux(double( or( ff<value_min, shift(ff,1,1)<value_min) ))=0.e0;
      flagno= and(ff>=value_min, shift(ff,1,1)>=value_min);
      udx=abs(xmassflux).*flagno;
      u2dt=DIVISION(dti2.*xmassflux.^2 .*flagno,(aru_3d.*AXB(dt_3d)));
      xmassflux=(udx-u2dt).*(DXB(ff).*AXB(dx_3d)./(2*AXB(ff)+epsilon)).*sw;
      xmassflux(double( lt(abs(udx),abs(u2dt)) ))=0.e0;
     
      ymassflux(double(or( ff<value_min, shift(ff,1,2)<value_min)))=0.e0;
      flagno= and(ff>=value_min, shift(ff,1,2)>=value_min);
      vdy=abs(ymassflux).*flagno;
      v2dt=dti2.*ymassflux.^2 ./(arv_3d.*AYB(dt_3d)).*flagno;
      ymassflux=(vdy-v2dt).*(DYB(ff).*AYB(dy_3d)./(2*AYB(ff)+epsilon)).*sw;
      ymassflux(double (lt(abs(vdy),abs(v2dt))))=0.e0;
      
      zwflux(double (or( ff<value_min, shift(ff,1,3)<value_min)))=0.e0;
      flagno= and(ff>=value_min, shift(ff,1,3)>=value_min);
      wdz=abs(zwflux) .* flagno;
      w2dt=DIVISION( dti2.*zwflux.^2 .* flagno , circshift(dzz_3d,1,3).*dt_3d );
      zwflux=(wdz-w2dt).*( -DZB(ff).* AZB(dz_3d) ./(2.*AZB(ff)+epsilon) ).*sw;
      zwflux(double (lt(abs(wdz),abs(w2dt))) )=0.e0;

end
