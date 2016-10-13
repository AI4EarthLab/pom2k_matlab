function F=DYF2_YZ(Y)
load('operator.mat');
F=OP_DYF2_YZ * Y * OP_R_YZ;
end
