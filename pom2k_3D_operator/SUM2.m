function F=SUM2(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(mx,kz);
for j=1:ny
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP_SUMZ2;
end
% for k=1:kz
%     F(:,:,k)=F(:,:,kz);
% end
end
