function F=DXF_XY(X)
load('operator.mat');
F=OP_DXF_XY * X * OP_R_XY;
end
