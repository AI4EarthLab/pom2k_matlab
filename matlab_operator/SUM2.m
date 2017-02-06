function F=SUM2(X)
global OP
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(mx,kz);
for j=1:ny
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP.OP_SUMZ2;
end
% for k=1:kz
%     F(:,:,k)=F(:,:,kz);
% end
end
