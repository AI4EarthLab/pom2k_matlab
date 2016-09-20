clear all;

load('test.mat');
load('operator.mat');

drhox=zeros(im,jm,kbm1+1);

%drhox(:,:,1)= grav*(-zz(1)) * AXB_XY(dt) .* DXB_XY(rho(:,:,1));


drhox1=drhox;
for k=2:kbm1
    for j=2:jmm1
        for i=2:imm1
            drhox(i,j,k)=drhox(i,j,k-1) ...
                         +grav*.25e0*(zz(k-1)-zz(k))*(dt(i,j)+dt(i-1,j))           ...
                          *(rho(i,j,k)+rho(i,j,k-1)-rho(i-1,j,k)-rho(i-1,j,k-1))    ...
                         +grav*.25e0*(zz(k-1)+zz(k))*(dt(i,j)-dt(i-1,j))           ...
                          *(rho(i,j,k)+rho(i-1,j,k)-rho(i,j,k-1)-rho(i-1,j,k-1));
        end
    end
end


rho_tmp=zeros(im,kb);
for j=2:jmm1
    % get a 2D array from a 3D array. 
    rho_tmp(:,:)=rho(:,j,:); 
    drhox1(:,j,:) = SUMZ_XZ(-grav * DZB_XZ(repmat(zz,im,1)) .* AXB_XZ(repmat(dt(:,j),1,kb)) .* DXB_XZ( AZB_XZ(rho_tmp) ) ...
                     +grav * AZB_XZ(repmat(zz,im,1)) .* DXB_XZ(repmat(dt(:,j),1,kb)) .* DZB_XZ( AXB_XZ(rho_tmp) ));
end
