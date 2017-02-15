function F=SUM2(X)
global OP
dim=size(X);
F=zeros(dim(1),dim(2),dim(3));
tmp=zeros(dim(1),dim(3));
for j=1:dim(2)
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP.OP_SUMZ2;
end
% for k=1:kz
%     F(:,:,k)=F(:,:,kz);
% end
end
