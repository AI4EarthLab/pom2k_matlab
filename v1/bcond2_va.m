function [vaf] = bcond2_va(vaf,vabn,vabs,rfn,rfs,el,eln,els,ramp)
    global im imm1 jm jmm1 dvm grav h;
    tmp1 = zeros(im, jm); 	tmp2 = zeros(im, jm);
    tmp1(:, jmm1) = vabn;   tmp1(:, 2) = vabs;
	tmp2(:, jmm1) = eln;    tmp2(:, 2) = els;

	tmp = ramp * (tmp1 + rfn * sqrt(grav./h) .* (el - tmp2));
	vaf(2:imm1, jm) = tmp(2:imm1, jmm1);
	tmp = ramp * (tmp1 - rfs * sqrt(grav./h) .* (el - tmp2));
	vaf(2:imm1, 2) = tmp(2:imm1, 2); vaf(2:imm1, 1) = tmp(2:imm1, 2);
    vaf(im, 2:jmm1) = 0.e0; vaf(1, 2:jmm1) = 0.e0;
    vaf = vaf .* dvm;
    return