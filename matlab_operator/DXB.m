function F=DXB(X)
%load('operator.mat');
global OP
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=OP.OP_DXB*X(:,:,k);
end

end
