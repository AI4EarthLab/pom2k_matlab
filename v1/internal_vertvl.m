function [w]=internal_vertvl (dt_3d,u,v,vfluxb,vfluxf,etf,etb)

global  im kb dz_3d dti2 fsm_3d;
w=( repmat(0.5e0*(vfluxb+vfluxf),1,1,kb)+SUM2(dz_3d .* ( ( DXF(AXB(dt_3d).* u)+DYF(AYB(dt_3d).* v) ) ...
    + repmat((etf-etb)/dti2,1,1,kb))) ).*fsm_3d;  %by hx
w(1,:,:) = 0.0; w(im,:,:) = 0.0;
end