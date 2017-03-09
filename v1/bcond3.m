function [uf, vf] = bcond3(uf, u, vf, v)
  
    global im imm1 jm jmm1 kb kbm1 hmax dum_3d dvm_3d h;

    tmph = h; tmph(2, :) = tmph(1, :); tmph = repmat(tmph, 1, 1, kb);
	tmpu = u; tmpu(im, :, :) = tmpu(imm1, :, :); tmpu(2, :, :) = tmpu(3, :, :);

	tmp = sqrt(tmph/hmax) .* AYF(AYB(tmpu)) + (1.0 - sqrt(tmph/hmax)) .* AYF(AYB(u));
	uf(im, 2:jmm1, 1:kbm1) = tmp(im, 2:jmm1, 1:kbm1);
    uf(2, 2:jmm1, 1:kbm1) = tmp(2, 2:jmm1, 1:kbm1);
	uf(1, 2:jmm1, 1:kbm1) = tmp(2, 2:jmm1, 1:kbm1);
    uf(2:imm1, jm, 1:kbm1) = 0.e0;
    uf(2:imm1, 1, 1:kbm1) = 0.e0;
    
    tmph = h; tmph(:, 2) = tmph(:, 1); tmph = repmat(tmph, 1, 1, kb);
    tmpv = v; tmpv(:, jm, :) = tmpv(:, jmm1, :); tmpv(:, 2, :) = tmpv(:, 3, :);

    tmp = sqrt(tmph/hmax) .* AXF(AXB(tmpv)) + (1.0 - sqrt(tmph/hmax)) .* AXF(AXB(v));
    vf(2:imm1, jm, 1:kbm1) = tmp(2:imm1, jm, 1:kbm1);
    vf(2:imm1, 2, 1:kbm1) = tmp(2:imm1, 2, 1:kbm1);
    vf(2:imm1, 1, 1:kbm1) = tmp(2:imm1, 2, 1:kbm1);
    vf(im, 2:jmm1, 1:kbm1) = 0.e0;
    vf(1, 2:jmm1, 1:kbm1) = 0.e0;

    uf = uf .* dum_3d;
    vf = vf .* dvm_3d;	
    
	return
end
