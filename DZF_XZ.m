function F=DZF_XZ(Z)
load('operator.mat');
F=OP_L_XZ * Z * OP_DZF_XZ;
end
