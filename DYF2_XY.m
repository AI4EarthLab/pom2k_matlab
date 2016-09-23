function F=DYF2_XY(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DYF2_XY;
end
