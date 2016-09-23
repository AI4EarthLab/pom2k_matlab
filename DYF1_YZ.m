function F=DYF1_YZ(Y)
load('operator.mat');
F=OP_DYF1_YZ * Y * OP_R_YZ;
end
