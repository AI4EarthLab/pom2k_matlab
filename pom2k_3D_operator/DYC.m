function F=DYC(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=X(:,:,k)*(OP_AYF * OP_DYB);
end

end