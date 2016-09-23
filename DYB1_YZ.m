function F=DYB1_YZ(Y)
load('operator.mat');
F=OP_DYB1_YZ * Y *OP_R_YZ;
end
