function F=AXB(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=OP_AXB1_XY*X(:,:,k);
end

end
