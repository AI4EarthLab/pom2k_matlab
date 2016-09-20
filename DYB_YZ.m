function F=DYB_YZ(Y)
load('operator.mat');
F=OP_DYB_YZ * Y *OP_R_YZ;
end
