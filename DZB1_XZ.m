function F=DZB1_XZ(Z)
load('operator.mat');
F=OP_L_XZ * Z * OP_DZB1_XZ;
end
