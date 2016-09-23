function F=DYB1_XY(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DYB1_XY;
end
