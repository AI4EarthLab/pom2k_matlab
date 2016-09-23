function F=DXF1_XZ(X)
load('operator.mat');
F=OP_DXF1_XZ * X * OP_R_XZ;
end
