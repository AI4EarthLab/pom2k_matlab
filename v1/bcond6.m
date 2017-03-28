function [uf, vf] = bcond6(uf, vf, u, v, q2, q2l)
    global im jm dti dx_3d dy_3d small fsm_3d fsm;
    tmpu = u;               tmpu(1, :, :) = tmpu(2, :, :);
	tmpv = v;               tmpv(:, 1, :) = tmpv(:, 2, :);
    ewflag=any(fsm,2);      snflag=any(fsm,1);

 %------------------EAST
 if(ewflag(im))   
    tmp1 = q2 - dti * ((.5e0 * (u - abs(u))) .* (small - q2) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(q2));
    uf(im, :, :) = tmp1(im, :, :);
    
	tmp1 = q2l - dti * ((.5e0 * (u - abs(u))) .* (small - q2l) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(q2l));
    vf(im, :, :) = tmp1(im, :, :);
 end    

 %------------------WEST 
 if(ewflag(1))
    tmp1 = q2 - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (q2 - small) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(q2));
	uf(1, :, :) = tmp1(1, :, :);

    tmp1 = q2l - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (q2l - small) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(q2l));
	vf(1, :, :) = tmp1(1, :, :);        
 end

 %------------------NORTH
 if(snflag(jm)) 
    tmp1 = q2 - dti * ((.5e0 * (v - abs(v))) .* (small - q2) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(q2));
	uf(:, jm, :) = tmp1(:, jm, :);
    
    tmp1 = q2l - dti * ((.5e0 * (v - abs(v))) .* (small - q2l) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(q2l));
	vf(:, jm, :) = tmp1(:, jm, :);        
 end

%------------------SOUTH 
 if(snflag(1)) 
	tmp1 = q2 - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (q2 - small) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(q2));
	uf(:, 1, :) = tmp1(:, 1, :);
 
	tmp1 = q2l - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (q2l - small) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(q2l));
	vf(:, 1, :) = tmp1(:, 1, :);    
 end

	uf = uf .* fsm_3d + 1.e-10;
	vf = vf .* fsm_3d + 1.e-10;
	return
