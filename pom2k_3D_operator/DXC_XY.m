function F=DXC_XY(X)
load('operator.mat');
F=OP_DXB2_XY * OP_AXF1_XY * X;
end