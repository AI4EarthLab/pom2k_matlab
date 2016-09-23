function F=DXB2_XZ(X)
load('operator.mat');
F=OP_DXB2_XZ * X *OP_R_XZ;
end
