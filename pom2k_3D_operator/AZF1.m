function F=AZF1(X)
%load('operator.mat');
global OP
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(ny,kz);
for i=1:mx
    tmp(:,:)=X(i,:,:);
    F(i,:,:)=tmp*OP.OP_AZF;
end

end