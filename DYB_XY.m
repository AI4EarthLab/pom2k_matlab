function F=DYB_XY(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DYB_XY;
end
