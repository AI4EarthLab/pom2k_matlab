function [uaf, vaf] = bcond2(uaf,uabe,uabw,rfe,rfw,vaf,vabn,vabs,rfn,rfs,el,ele,eln,elw,els,ramp)
    global im imm1 jm jmm1 dum dvm grav h;

    tmp1 = zeros(im, jm); tmp1(imm1, :) = uabe;     tmp1(2, :) = uabw;
	tmp2 = zeros(im, jm); tmp2(imm1, :) = ele;      tmp2(2, :) = elw;
 
	tmp = ramp * (tmp1 + rfe * sqrt(grav./h) .* (el - tmp2));
	uaf(im, 2:jmm1) = tmp(imm1, 2:jmm1);
	tmp = ramp * (tmp1 - rfw * sqrt(grav./h) .* (el - tmp2));
	uaf(2, 2:jmm1) = tmp(2, 2:jmm1); uaf(1, 2:jmm1) = tmp(2, 2:jmm1);
	uaf(2:imm1, jm) = 0.e0; uaf(2:imm1, 1) = 0.e0;
    uaf = uaf .* dum;

    tmp1(:, jmm1) = vabn; tmp1(:, 2) = vabs;
	tmp2(:, jmm1) = eln; tmp2(:, 2) = els;

	tmp = ramp * (tmp1 + rfn * sqrt(grav./h) .* (el - tmp2));
	vaf(2:imm1, jm) = tmp(2:imm1, jmm1);
	tmp = ramp * (tmp1 - rfs * sqrt(grav./h) .* (el - tmp2));
	vaf(2:imm1, 2) = tmp(2:imm1, 2); vaf(2:imm1, 1) = tmp(2:imm1, 2);
    vaf(im, 2:jmm1) = 0.e0; vaf(1, 2:jmm1) = 0.e0;
    vaf = vaf .* dvm;
    return
end