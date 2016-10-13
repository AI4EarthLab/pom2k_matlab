function F=DZC_XZ(Z)
load('operator.mat');
F=Z*(OP_AZF1_XZ * OP_DZB2_XZ);
end