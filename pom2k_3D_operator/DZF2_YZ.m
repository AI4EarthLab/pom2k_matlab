function F=DZF2_YZ(Z)
load('operator.mat');
F=OP_L_YZ * Z * OP_DZF2_YZ;
end
