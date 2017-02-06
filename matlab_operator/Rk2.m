function F=Rk2(X)
%load('operator.mat');
global OP
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(mx,kz);
for j=1:ny
    tmp(:,:)=X(:,j,:);
    F(:,j,:)=tmp*OP.OP_Rk;
end

end