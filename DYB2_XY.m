function F=DYB2_XY(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DYB2_XY;
end
