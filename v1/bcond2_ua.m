function [uaf] = bcond2_ua(uaf,uabe,uabw,rfe,rfw,el,ele,elw,ramp)
    global im imm1 jm jmm1 dum grav h;

    tmp1 = zeros(im, jm); tmp1(imm1, :) = uabe;     tmp1(2, :) = uabw;
	tmp2 = zeros(im, jm); tmp2(imm1, :) = ele;      tmp2(2, :) = elw;
 
	tmp = ramp * (tmp1 + rfe * sqrt(grav./h) .* (el - tmp2));
	uaf(im, 2:jmm1) = tmp(imm1, 2:jmm1);
	tmp = ramp * (tmp1 - rfw * sqrt(grav./h) .* (el - tmp2));
	uaf(2, 2:jmm1) = tmp(2, 2:jmm1); uaf(1, 2:jmm1) = tmp(2, 2:jmm1);
	uaf(2:imm1, jm) = 0.e0; uaf(2:imm1, 1) = 0.e0;
    uaf = uaf .* dum;
return