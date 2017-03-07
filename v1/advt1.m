function [ff]=advt1(fb,f,dt_3d,u,v,aam,w,etb,etf)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Integrates conservative scalar equations.           *
% *                                                                    *
% *                This is centred scheme, as originally provide in    *
% *                POM (previously called advt).                       *
% *                                                                    *
% **********************************************************************
global kb dum_3d dvm_3d dti2 tprni h_3d;

ff= ( (h_3d+repmat(etb,1,1,kb)).*fb - dti2 .* (DXF( AXB(dt_3d).*AXB(f).*u-AXB(aam).*AXB(h_3d)*tprni.*DXB(fb).*dum_3d ) ...
    + DYF( AYB(dt_3d).*AYB(f).*v-AYB(aam).*AYB(h_3d)*tprni.*DYB(fb).*dvm_3d )-DZF( AZB(f).*w )) ) ./((h_3d+repmat(etf,1,1,kb))) ; 
end