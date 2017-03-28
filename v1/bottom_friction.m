function [cbc] = bottom_friction(kappa,zz,h,z0b,cbcmin,cbcmax)
global kbm1;
cbc=(kappa./log((1.0+zz(kbm1))*h/z0b)).^2;
cbc=max(cbcmin,cbc);        cbc=min(cbcmax,cbc);

end