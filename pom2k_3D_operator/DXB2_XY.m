function F=DXB2_XY(X)
load('operator.mat');
F=OP_DXB2_XY * X * OP_R_XY;
end
