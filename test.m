clear all;

load('test.mat');
load('operator.mat');

drhoy=zeros(im,jm,kbm1+1);
A=zeros(im,jm,kbm1+1);
B=A;
%drhox(:,:,1)= grav*(-zz(1)) * AXB_XY(dt) .* DXB_XY(rho(:,:,1));

%
drhoy1=drhoy;

drhoy(:,:,1)= grav*(-zz(1)) * AYB_XY(dt) .* DYB_XY(rho(:,:,1));
%
for k=2:kbm1
    for j=2:jmm1
        for i=2:imm1
            A(i,j,k)=0.5*(dt(i,j)+dt(i,j-1));
            drhoy(i,j,k)=drhoy(i,j,k-1)+...
                          +grav*0.25*(zz(k-1)-zz(k))*(dt(i,j)+dt(i,j-1))            ...
                           *(rho(i,j,k)-rho(i,j-1,k)+rho(i,j,k-1)-rho(i,j-1,k-1))     ...
                          +grav*0.25*(zz(k-1)+zz(k))*(dt(i,j)-dt(i,j-1))            ...
                           *(rho(i,j,k)+rho(i,j-1,k)-rho(i,j,k-1)-rho(i,j-1,k-1));
        end
    end
end

rho_tmp=zeros(jm,kb);
for i=2:imm1
    % get a 2D array from a 3D array. 
    rho_tmp(:,:)=rho(i,:,:); 
    tmp=repmat(drhoy(i,:,1)',1,kb) ;
    tmp(:,kb)=0;
    B(i,:,:)=AYB_YZ(repmat(dt(i,:)',1,kb));
    drhoy1(i,:,:) = tmp+ SUMZ_YZ(-grav * DZB_YZ(repmat(zz,jm,1)) .* AYB_YZ(repmat(dt(i,:)',1,kb)) .* DYB_YZ( AZB_YZ(rho_tmp) ) ...
                                   +grav * AZB_YZ(repmat(zz,jm,1)) .* DYB_YZ(repmat(dt(i,:)',1,kb)) .* DZB_YZ( AYB_YZ(rho_tmp) )); 
end
