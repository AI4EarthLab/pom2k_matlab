function F=DZB2_XZ(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DZB2_XZ;
end
