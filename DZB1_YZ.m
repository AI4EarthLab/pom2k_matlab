function F=DZB1_YZ(Z)
load('operator.mat');
F=OP_L_YZ * Z * OP_DZB1_YZ;
end
