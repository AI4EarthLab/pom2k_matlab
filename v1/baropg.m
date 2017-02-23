 function [rho,drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
% *                                                                    *
% **********************************************************************
global zz_3d grav dum_3d dvm_3d dz_3d

rho = rho - rmean;
% Calculate x-component of baroclinic pressure gradient:
% % % drhox(:,:,1)= -grav .* zz_3d(:,:,1) .* AXB(dt_3d(:,:,1)) .* DXB(rho(:,:,1));
% % % drhox= SUM1(drhox-grav *DZB(zz_3d) .* AXB(dt_3d) .* DXB(AZB(rho)) + grav * AZB(zz_3d) .* DXB(dt_3d) .* DZB(AXB(rho)));
% % % drhox = ramp*drhox .* AXB(dt_3d) .* dum_3d .* AXB(dy_3d);
drhox=ramp * grav * AXB(dt_3d) .* SUMZ(-DZB(zz_3d).* AXB(dt_3d) .* DXB(AZB(rho)) + AZB(zz_3d) .* DXB(dt_3d).* DZB(AXB(rho)).*AZB(dz_3d)).* dum_3d;

%     Calculate y-component of baroclinic pressure gradient:
% drhoy(:,:,1)= -grav .* zz_3d(:,:,1) .* AYB(dt_3d(:,:,1)) .* DYB(rho(:,:,1));
% drhoy= SUM1(drhoy-grav * DZB(zz_3d) .* AYB(dt_3d) .* DYB(AZB(rho)) + grav * AZB(zz_3d) .* DYB(dt_3d) .* DZB(AYB(rho)));
% drhoy = ramp*drhoy .* AYB(dt_3d) .* dvm_3d .* AYB(dx_3d);

drhoy=ramp * grav .* AYB(dt_3d) .* SUMZ(-DZB(zz_3d) .* AYB(dt_3d) .* DYB(AZB(rho))+ AZB(zz_3d) .* DYB(dt_3d) .* DZB(AYB(rho)).*AZB(dz_3d)).* dvm_3d;
%
rho = rho + rmean;