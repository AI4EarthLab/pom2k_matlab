function F=DXF1_XY(X)
load('operator.mat');
F=OP_DXF1_XY * X * OP_R_XY;
end
