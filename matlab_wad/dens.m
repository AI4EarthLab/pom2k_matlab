function [si,ti,rhoo]=dens(si,ti,rhoo,pdens, ...
    im,jm,kbm1,tbias,sbias,grav,rhoref,zz,h,fsm) 
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates (density-1000.)/rhoref.                  *
% *                                                                    *
% *                (see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech.,  *
% *                609-611.)                                           *
% *                                                                    *
% *                ti is potential temperature                         *
% *                                                                    *
% *                If using 32 bit precision, it is recommended that   *
% *                cr,p,rhor,sr,tr,tr2,tr3 and tr4 be made double      *
% *                precision, and the "e"s in the constants be changed *
% *                to "d"s.                                            *
% *                                                                    *
% * NOTE: if pressure is not used in dens, buoyancy term (boygr)       *
% *       in profq must be changed (see note in profq)                 *
% *                                                                    *
% **********************************************************************


for k=1:kbm1
	for j=1:jm
    	for i=1:im
        	tr=ti(i,j,k)+tbias;
            sr=si(i,j,k)+sbias;
            tr2=tr*tr;
            tr3=tr2*tr;
            tr4=tr3*tr;
%     Approximate pressure in units of bars:
% WAD Zhaowei 2017/11/9
% lyo:!wad:Keep "p" defined as in original pom w/o hhi (i.e. "h" defined
%          wrt MSL), but set p=0 for cells where input "h" indicates
%          "dry" (i.e. where h-hhi < 0); i.e. use "pdens":
% 
            p=pdens(i,j,k);
%           p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5;
            rhor=-0.157406e0+6.793952e-2*tr-9.095290e-3*tr2+1.001685e-4*tr3...
				 -1.120083e-6*tr4+6.536332e-9*tr4*tr;
 
            rhor=rhor+(0.824493e0-4.0899e-3*tr+7.6438e-5*tr2-8.2467e-7*tr3...
                  +5.3875e-9*tr4)*sr+(-5.72466e-3+1.0227e-4*tr...
                  -1.6546e-6*tr2)*abs(sr)^1.5+4.8314e-4*sr*sr;
 
            cr=1449.1e0+.0821e0*p+4.55e0*tr-.045e0*tr2+1.34e0*(sr-35.e0);
            rhor=rhor+1.e5*p/(cr*cr)*(1.e0-2.e0*p/(cr*cr));
 
            rhoo(i,j,k)=rhor/rhoref*fsm(i,j);
		end
	end
end

