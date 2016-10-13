function F=DYC_XY(Y)
load('operator.mat');
F=Y*(OP_AYF1_XY * OP_DYB2_XY);
end