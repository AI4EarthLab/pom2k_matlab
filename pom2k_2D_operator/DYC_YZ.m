function F=DYC_YZ(Z)
load('operator.mat');
F=Z*(OP_AZF1_YZ * OP_DZB2_YZ);
end