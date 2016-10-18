function [rho,drhox,drhoy] = new_baropg(rho, rmean, dt, ramp)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('para.mat');

rho = rho - rmean;
drhox=zeros(im,jm,kb); drhoy=zeros(im,jm,kb);

dt_3d=repmat(dt,1,1,kb);
%     Calculate x-component of baroclinic pressure gradient:
drhox(:,:,1)= -grav .* zz_3d(:,:,1) .* AXB1(dt_3d(:,:,1)) .* DXB2(rho(:,:,1));
drhox= SUM1(drhox-grav *DZB2(zz_3d) .* AXB1(dt_3d) .* DXB2(AZB1(rho)) + grav * AZB1(zz_3d) .* DXB2(dt_3d) .* DZB2(AXB1(rho)));
drhox = ramp*drhox .* AXB1(dt_3d) .* dum_3d .* AXB1(dy_3d);

%     Calculate y-component of baroclinic pressure gradient:
drhoy(:,:,1)= -grav .* zz_3d(:,:,1) .* AYB1(dt_3d(:,:,1)) .* DYB2(rho(:,:,1));
drhoy= SUM1(drhoy-grav * DZB2(zz_3d) .* AYB1(dt_3d) .* DYB2(AZB1(rho)) + grav * AZB1(zz_3d) .* DYB2(dt_3d) .* DZB2(AYB1(rho)));
drhoy = ramp*drhoy .* AYB1(dt_3d) .* dvm_3d .* AYB1(dx_3d);

rho = rho + rmean;