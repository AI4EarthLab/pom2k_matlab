function F=DXC_XZ(X)
load('operator.mat');
F=OP_DXB2_XZ * OP_AXF1_XZ * X;
end