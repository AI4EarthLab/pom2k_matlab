function F=DXF2(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
    F(:,:,k)=OP_DXF2_XY*X(:,:,k)*OP_R_XY;
end

end
