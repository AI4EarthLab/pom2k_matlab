function F=DXF(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
%    F(:,:,k)=OP_DXF2_XY*X(:,:,k)*OP_R_XY;
    F(:,:,k)=OP_DXF1_XY*X(:,:,k);    
end

end
