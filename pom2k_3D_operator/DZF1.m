function F=DZF1(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(mx,kz);
for j=1:ny
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP_DZF2_XZ;
end

end
