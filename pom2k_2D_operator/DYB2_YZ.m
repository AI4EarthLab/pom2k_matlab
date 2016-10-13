function F=DYB2_YZ(Y)
load('operator.mat');
F=OP_DYB2_YZ * Y *OP_R_YZ;
end
