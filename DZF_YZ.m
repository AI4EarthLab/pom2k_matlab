function F=DZF_YZ(Z)
load('operator.mat');
F=OP_L_YZ * Z * OP_DZF_YZ;
end
