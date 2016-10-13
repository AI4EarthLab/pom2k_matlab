function F=DZB2_YZ(Z)
load('operator.mat');
F=OP_L_YZ * Z * OP_DZB2_YZ;
end
