function F=DYF_XY(Y)
load('operator.mat');
F=OP_L_XY * Y * OP_DYF_XY;
end
