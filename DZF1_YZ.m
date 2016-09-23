function F=DZF1_YZ(Z)
load('operator.mat');
F=OP_L_YZ * Z * OP_DZF1_YZ;
end
