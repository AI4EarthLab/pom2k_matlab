function F=DXB_XY(X)
load('operator.mat');
F=OP_DXB_XY * X * OP_R_XY;
end
