function F=AXF2(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=OP_AXF2_XY*X(:,:,k);
end

end
