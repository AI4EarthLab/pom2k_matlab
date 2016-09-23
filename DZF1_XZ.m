function F=DZF1_XZ(Z)
load('operator.mat');
F=OP_L_XZ * Z * OP_DZF1_XZ;
end
