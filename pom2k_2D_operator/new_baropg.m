function [rho,drhox,drhoy] = new_baropg(rho, rmean, dt, ramp)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
% *                                                                    *
% **********************************************************************
% 
load('grid.mat'); load('depth.mat'); load('para.mat'); load('masks.mat'); load('xyz.mat');

rho = rho - rmean;
drhox=zeros(im,jm,kb); drhoy=zeros(im,jm,kb);
%
%     Calculate x-component of baroclinic pressure gradient:
%
drhox(:,:,1)= grav*(-zz(1)) * AXB1_XY(dt) .* DXB2_XY(rho(:,:,1));
                    
rho_tmp=zeros(im,kb);
for j=1:jm
    % get a 2D array from a 3D array. 
    rho_tmp(:,:)=rho(:,j,:); 
    tmp=repmat(drhox(:,j,1),1,kb) ;
    tmp(:,kb)=0;
    drhox(:,j,:) =  tmp + SUMZ1(-grav * DZB2_XZ(repmat(zz,im,1)) .* AXB1_XZ(repmat(dt(:,j),1,kb)) .* DXB2_XZ( AZB1_XZ(rho_tmp) ) ...
                                  +grav * AZB1_XZ(repmat(zz,im,1)) .* DXB2_XZ(repmat(dt(:,j),1,kb)) .* DZB2_XZ( AXB1_XZ(rho_tmp) )); 
end

%save('test.mat','im','jm','kb','kbm1','jmm1','imm1','grav','zz','dt','rho');

for k=1:kbm1
    drhox(:,:,k) = drhox(:,:,k) .* AXB1_XY(dt) .* dum .* AXB1_XY(dy);
end
 
%
%     Calculate y-component of baroclinic pressure gradient:
%
drhoy(:,:,1)= grav*(-zz(1)) * AYB1_XY(dt) .* DYB2_XY(rho(:,:,1));
%
rho_tmp=zeros(jm,kb);
for i=2:imm1
    % get a 2D array from a 3D array. 
    rho_tmp(:,:)=rho(i,:,:); 
    tmp=repmat(drhoy(i,:,1)',1,kb) ;
    tmp(:,kb)=0;
    drhoy(i,:,:) =  tmp + SUMZ1(-grav * DZB2_YZ(repmat(zz,jm,1)) .* AYB1_YZ(repmat(dt(i,:)',1,kb)) .* DYB2_YZ( AZB1_YZ(rho_tmp) ) ...
                                   +grav * AZB1_YZ(repmat(zz,jm,1)) .* DYB2_YZ(repmat(dt(i,:)',1,kb)) .* DZB2_YZ( AYB1_YZ(rho_tmp) )); 
end

for k=1:kbm1
    drhoy(:,:,k) = drhoy(:,:,k) .* AYB1_XY(dt) .* dvm .* AYB1_XY(dx);
end

drhox=ramp*drhox;
drhoy=ramp*drhoy;

rho = rho + rmean;