function F=DXF2_XY(X)
load('operator.mat');
F=OP_DXF2_XY * X * OP_R_XY;
end
