function F=DXF2_XZ(X)
load('operator.mat');
F=OP_DXF2_XZ * X * OP_R_XZ;
end
