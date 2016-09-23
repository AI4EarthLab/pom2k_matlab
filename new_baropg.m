function [rho,drhox,drhoy] = new_baropg(rho, rmean, dt, ramp)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
% *                                                                    *
% **********************************************************************
% 
load('grid.mat');
load('depth.mat');
load('para.mat');
load('masks.mat');
load('operator.mat');
load('xyz.mat');

rho = rho - rmean;
drhox=zeros(im,jm,kb); drhoy=zeros(im,jm,kb);
%
%     Calculate x-component of baroclinic pressure gradient:
%
drhox(:,:,1)= grav*(-zz(1)) * AXB_XY(dt) .* DXB_XY(rho(:,:,1));
                    
rho_tmp=zeros(im,kb);
for j=1:jm
    % get a 2D array from a 3D array. 
    rho_tmp(:,:)=rho(:,j,:); 
    tmp=repmat(drhox(:,j,1),1,kb) ;
    tmp(:,kb)=0;
    drhox(:,j,:) =  tmp + SUMZ_XZ(-grav * DZB_XZ(repmat(zz,im,1)) .* AXB_XZ(repmat(dt(:,j),1,kb)) .* DXB_XZ( AZB_XZ(rho_tmp) ) ...
                                  +grav * AZB_XZ(repmat(zz,im,1)) .* DXB_XZ(repmat(dt(:,j),1,kb)) .* DZB_XZ( AXB_XZ(rho_tmp) )); 
end

%save('test.mat','im','jm','kb','kbm1','jmm1','imm1','grav','zz','dt','rho');

for k=1:kbm1
    drhox(:,:,k) = drhox(:,:,k) .* AXB_XY(dt) .* dum .* AXB_XY(dy);
end
 
%
%     Calculate y-component of baroclinic pressure gradient:
%
drhoy(:,:,1)= grav*(-zz(1)) * AYB_XY(dt) .* DYB_XY(rho(:,:,1));
%
rho_tmp=zeros(jm,kb);
for i=2:imm1
    % get a 2D array from a 3D array. 
    rho_tmp(:,:)=rho(i,:,:); 
    tmp=repmat(drhoy(i,:,1)',1,kb) ;
    tmp(:,kb)=0;
    drhoy(i,:,:) =  tmp + SUMZ_YZ(-grav * DZB_YZ(repmat(zz,jm,1)) .* AYB_YZ(repmat(dt(i,:)',1,kb)) .* DYB_YZ( AZB_YZ(rho_tmp) ) ...
                                   +grav * AZB_YZ(repmat(zz,jm,1)) .* DYB_YZ(repmat(dt(i,:)',1,kb)) .* DZB_YZ( AYB_YZ(rho_tmp) )); 
end

for k=1:kbm1
    drhoy(:,:,k) = drhoy(:,:,k) .* AYB_XY(dt) .* dvm .* AYB_XY(dx);
end

drhox=ramp*drhox;
drhoy=ramp*drhoy;

rho = rho + rmean;