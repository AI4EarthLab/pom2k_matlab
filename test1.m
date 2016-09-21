clear all;

load('test.mat');
load('operator.mat');
load('AB.mat');

drhox=zeros(im,jm,kbm1+1);
%drhox(:,:,1)= grav*(-zz(1)) * AXB_XY(dt) .* DXB_XY(rho(:,:,1));


drhox1=drhox;
for k=2:kbm1
    for j=2:jmm1
        for i=2:imm1
            drhox(i,j,k)=drhox(i,j,k-1)+A(i,j,k);
        end
    end
end

R=triu(ones(kb,kb)) ;

OP_SUMZ_XZ= R;

tmp=zeros(im,kb);
for j=2:jmm1
    tmp(:,:)=A(:,j,:);
    drhox1(:,j,:) = tmp(:,:)*OP_SUMZ_XZ;
end