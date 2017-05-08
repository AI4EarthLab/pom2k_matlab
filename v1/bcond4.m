function [tf] = bcond4(tf,u,v,w,t,dt,tbe,tbw,tbn,tbs)
    global im imm1 jm jmm1 kb kbm1 kbm2 dti fsm dx_3d dy_3d fsm_3d dz_3d dzz_3d;

	tmptb = zeros(im, jm, kb); tmptb(im, :, :) = tbe; tmptb(1, :, :) = tbw; tmptb(:, jm, :) = tbn; tmptb(:, 1, :) = tbs;
    tmpu = u; tmpu(imm1, :, :) = tmpu(im, :, :);    tmpu(1, :, :) = tmpu(2, :, :);
	tmpv = v; tmpv(:, jmm1, :) = tmpv(:, jm, :);    tmpv(:, 1, :) = tmpv(:, 2, :);
	dt_3d = repmat(dt, 1, 1, kb)               ;    ewflag=any(fsm,2);                             snflag=any(fsm,1);                 
%------------------EAST
 if(ewflag(im))
    tmp1 = t - dti * ((.5e0 * (u - abs(u))) .* (tmptb - t) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(t));
    tf(im, :, 1:kbm1) = tmp1(im, :, 1:kbm1);    tmp1(imm1, :, :) = tmp1(im, :, :);
    
    tmp2 = tmp1 + dti * DIVISION(tmpu + abs(tmpu), 2 * tmpu) .* AZF(w)./dt_3d.*DZF(AZB(t)).*DIVISION(dz_3d,AZB(dzz_3d));    %add by hx
    tf(im, :, 2:kbm2) = tmp2(imm1, :, 2:kbm2);
 end   
%------------------WEST 
 if(ewflag(1))
    tmp1 = t - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (t - tmptb) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(t).*dx_3d./(AXF(dx_3d)) );  %add by hx
	tf(1, :, 1:kbm1) = tmp1(1, :, 1:kbm1);      tmp1(2, :, :) = tmp1(1, :, :);
    
	tmp2 = tmp1 + dti * DIVISION(u - abs(u), 2 * u) .* AZF(w) ./ dt_3d .* DZF(AZB(t)).*DIVISION(dz_3d,AZB(dzz_3d));         %add by hx
	tf(1, :, 2:kbm2) = tmp2(2, :, 2:kbm2);
 end   
 %------------------NORTH
 if(snflag(jm))
    tmp1 = t - dti * ((.5e0 * (v - abs(v))) .* (tmptb - t) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(t));
	tf(:, jm, 1:kbm1) = tmp1(:, jm, 1:kbm1);	tmp1(:, jmm1, :) = tmp1(:, jm, :);
    
	tmp2 = tmp1 + dti * DIVISION(tmpv + abs(tmpv), 2 * tmpv) .* AZF(w) ./dt_3d .*DZF(AZB(t)).*DIVISION(dz_3d,AZB(dzz_3d));  %add by hx
	tf(:, jm, 2:kbm2) = tmp2(:, jmm1, 2:kbm2);
 end   
%------------------SOUTH 
 if(snflag(1))
	tmp1 = t - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (t - tmptb) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(t));
	tf(:, 1, 1:kbm1) = tmp1(:, 1, 1:kbm1);   	tmp1(:, 2, :) = tmp1(:, 1, :);
    
	tmp2 = tmp1 + dti * DIVISION(v - abs(v), 2 * v) .* AZF(w) ./ dt_3d .* DZF(AZB(t)).*DIVISION(dz_3d,AZB(dzz_3d));         %add by hx
	tf(:, 1, 2:kbm2) = tmp2(:, 2, 2:kbm2);
 end   
    
	tf = tf .* fsm_3d;
end