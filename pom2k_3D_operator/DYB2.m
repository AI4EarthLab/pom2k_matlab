function F=DYB2(X)
load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
for k=1:kz
%    F(:,:,k)=OP_L_XY *X(:,:,k)*OP_DYB2_XY;
    F(:,:,k)=X(:,:,k)*OP_DYB2_XY;    
end

end
