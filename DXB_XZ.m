function F=DXB_XZ(X)
load('operator.mat');
F=OP_DXB_XZ * X *OP_R_XZ;
end
