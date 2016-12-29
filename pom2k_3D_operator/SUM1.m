function F=SUM1(X)
%load('operator.mat');
global OP
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(mx,kz);
for j=1:ny
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP.OP_SUMZ1;
end
% for k=1:kz
%     F(:,:,k)=F(:,:,k);
% end
end
