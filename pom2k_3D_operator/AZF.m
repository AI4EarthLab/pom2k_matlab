function F=AZF(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(mx,kz);
for j=1:ny
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP_AZF1_XZ;
end

end
