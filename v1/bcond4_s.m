function [sf] = bcond4_s(sf, u, v, w, s, dt)
    global im imm1 jm jmm1 kb kbm1 kbm2 dti dx_3d dy_3d fsm_3d sbe sbw sbn sbs;

	tmpsb = zeros(im, jm, kb); tmpsb(im, :, :) = sbe; tmpsb(1, :, :) = sbw; tmpsb(:, jm, :) = sbn; tmpsb(:, 1, :) = sbs;
    tmpu = u; tmpu(imm1, :, :) = tmpu(im, :, :); tmpu(1, :, :) = tmpu(2, :, :);
	tmpv = v; tmpv(:, jmm1, :) = tmpv(:, jm, :); tmpv(:, 1, :) = tmpv(:, 2, :);
	dt_3d = repmat(dt, 1, 1, kb);
    
    tmp1 = s - dti * ((.5e0 * (u - abs(u))) .* (tmpsb - s) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(s));
    sf(im, :, 1:kbm1) = tmp1(im, :, 1:kbm1);
    sf(im, :, 1) = tmp1(im, :, 1); sf(im, :, kbm1) = tmp1(im, :, kbm1);
    tmp1(imm1, :, :) = tmp1(im, :, :);
    tmp2 = tmp1 - dti * DIVISION(tmpu + abs(tmpu), 2 * tmpu) .* AZF(w) ./ dt_3d	.* DZF(AZB(s));
    sf(im, :, 2:kbm2) = tmp2(imm1, :, 2:kbm2);

    tmp1 = s - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (s - tmpsb) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(s));
	sf(1, :, 1:kbm1) = tmp1(1, :, 1:kbm1);
    sf(1, :, 1) = tmp1(1, :, 1); sf(1, :, kbm1) = tmp1(1, :, kbm1);
	tmp1(2, :, :) = tmp1(1, :, :);
	tmp2 = tmp1 - dti * DIVISION(u - abs(u), 2 * u) .* AZF(w) ./ dt_3d .* DZF(AZB(s));
	sf(1, :, 2:kbm2) = tmp2(2, :, 2:kbm2);

    tmp1 = s - dti * ((.5e0 * (v - abs(v))) .* (tmpsb - s) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(s));
	sf(:, jm, 1:kbm1) = tmp1(:, jm, 1:kbm1);
    sf(:, jm, 1) = tmp1(:, jm, 1); sf(:, jm, kbm1) = tmp1(:, jm, kbm1);
	tmp1(:, jmm1, :) = tmp1(:, jm, :);
	tmp2 = tmp1 - dti * DIVISION(tmpv + abs(tmpv), 2 * tmpv) .* AZF(w) ./ dt_3d .* DZF(AZB(s));
	sf(:, jm, 2:kbm2) = tmp2(:, jmm1, 2:kbm2);

	tmp1 = s - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (s - tmpsb) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(s));
	sf(:, 1, 1:kbm1) = tmp1(:, 1, 1:kbm1);
    sf(:, 1, 1) = tmp1(:, 1, 1); sf(:, 1, kbm1) = tmp1(:, 1, kbm1);
	tmp1(:, 2, :) = tmp1(:, 1, :);
	tmp2 = tmp1 - dti * DIVISION(v - abs(v), 2 * v) .* AZF(w) ./ dt_3d .* DZF(AZB(s));
	sf(:, 1, 2:kbm2) = tmp2(:, 2, 2:kbm2);
	sf = sf .* fsm_3d;
	return
end