function F=DYB2(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=X(:,:,k)*OP_DYB;    
end

end
