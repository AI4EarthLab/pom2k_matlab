function F=DYF1_XY(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DYF1_XY;
end
