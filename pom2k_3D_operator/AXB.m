function F=AXB(X)
%load('operator.mat');
global OP
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);

for k=1:kz
    F(:,:,k)=OP.OP_AXB*X(:,:,k);
end

end
