function F=AXB1(X)
global OP
%load('operator.mat');
[mx,ny,kz]=size(X);
F=zeros(mx,ny,kz);
tmp=zeros(ny,kz);
for i=1:mx
    tmp(:,:)=X(i,:,:);
    F(i,:,:)=OP.OP_AXB1 * tmp;
end

end