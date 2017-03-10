function [sf,s,sb] = bcond4_s(sf,sb,u,v,w,s,dt)
    global im imm1 jm jmm1 kb zz kbm1 dti dx dx_3d dy dy_3d fsm_3d sbe sbw sbn sbs;

	tmpsb = zeros(im, jm, kb); tmpsb(im, :, :) = sbe; tmpsb(1, :, :) = sbw; tmpsb(:, jm, :) = sbn; tmpsb(:, 1, :) = sbs;
    tmpu = u; tmpu(imm1, :, :) = tmpu(im, :, :); tmpu(1, :, :) = tmpu(2, :, :);
	tmpv = v; tmpv(:, jmm1, :) = tmpv(:, jm, :); tmpv(:, 1, :) = tmpv(:, 2, :);

%-------------------------------------------------sf--------------------------------------------------------
    tmp1 = s - dti * ((.5e0 * (u - abs(u))) .* (tmpsb - s) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(s));
    sf(im, :, 1:kbm1) = tmp1(im, :, 1:kbm1);

    tmp1 = s - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (s - tmpsb) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(s));
	sf(1, :, 1:kbm1) = tmp1(1, :, 1:kbm1);

    tmp1 = s - dti * ((.5e0 * (v - abs(v))) .* (tmpsb - s) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(s));
	sf(:, jm, 1:kbm1) = tmp1(:, jm, 1:kbm1);

	tmp1 = s - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (s - tmpsb) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(s));
	sf(:, 1, 1:kbm1) = tmp1(:, 1, 1:kbm1);

    for k=1:kbm1
        for j=1:jm
            if(double(u(im, j, k) > 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(imm1,j));
                    sf(im,j,k)= s(im,j,k)-2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))*(s(im,j,k)-s(imm1,j,k))-wm*(s(imm1,j,k-1)-s(imm1,j,k+1));
                end
            end
            
            if(double(u(2, j, k) < 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(2,j));
                    sf(1,j,k)=s(1,j,k)-2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))*(s(2,j,k)-s(1,j,k))-wm*(s(2,j,k-1)-s(2,j,k+1));
                end
            end
        end
    end
    for k=1:kbm1
        for i=1:im           
            if(double(v(i, jm, k) > 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(i,jmm1));
                    sf(i,jm,k)=s(i,jm,k)-2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))*(s(i,jm,k)-s(i,jmm1,k))-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1));
                end
            end
            if(double(v(i, 2, k) < 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(i,2));
                    sf(i,1,k)=s(i,1,k)-2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))*(s(i,2,k)-s(i,1,k))-wm*(s(i,2,k-1)-s(i,2,k+1));
                end
            end
        end
    end

	sf = sf .* fsm_3d;

	return

