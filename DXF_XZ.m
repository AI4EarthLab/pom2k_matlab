function F=DXF_XZ(X)
load('operator.mat');
F=OP_DXF_XZ * X * OP_R_XZ;
end
