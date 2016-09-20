function F=DZB_XZ(Z)
load('operator.mat');
F=OP_L_XZ * Z * OP_DZB_XZ;
end
