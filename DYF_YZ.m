function F=DYF_YZ(Y)
load('operator.mat');
F=OP_DYF_YZ * Y * OP_R_YZ;
end
