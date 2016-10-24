function F=DXB(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=OP_DXB1_XY*X(:,:,k);
end

end
