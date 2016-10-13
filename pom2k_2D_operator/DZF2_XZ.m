function F=DZF2_XZ(Z)
load('operator.mat');
F=OP_L_XZ * Z * OP_DZF2_XZ;
end
