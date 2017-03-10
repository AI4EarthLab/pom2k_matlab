function [tf,t,tb] = bcond4_t(tf,tb,u,v,w,t,dt)
    global im imm1 jm jmm1 kb zz kbm1 dti dx dx_3d dy dy_3d fsm_3d tbe tbw tbn tbs;

	tmptb = zeros(im, jm, kb); tmptb(im, :, :) = tbe; tmptb(1, :, :) = tbw; tmptb(:, jm, :) = tbn; tmptb(:, 1, :) = tbs;
    tmpu = u; tmpu(imm1, :, :) = tmpu(im, :, :); tmpu(1, :, :) = tmpu(2, :, :);
	tmpv = v; tmpv(:, jmm1, :) = tmpv(:, jm, :); tmpv(:, 1, :) = tmpv(:, 2, :);

%-------------------------------------------------tf--------------------------------------------------------
    tmp1 = t - dti * ((.5e0 * (u - abs(u))) .* (tmptb - t) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(t));
    tf(im, :, 1:kbm1) = tmp1(im, :, 1:kbm1);
 
    tmp1 = t - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (t - tmptb) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(t));
	tf(1, :, 1:kbm1) = tmp1(1, :, 1:kbm1);
 
    tmp1 = t - dti * ((.5e0 * (v - abs(v))) .* (tmptb - t) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(t));
	tf(:, jm, 1:kbm1) = tmp1(:, jm, 1:kbm1);
    
	tmp1 = t - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (t - tmptb) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(t));
	tf(:, 1, 1:kbm1) = tmp1(:, 1, 1:kbm1);

    for k=1:kbm1
        for j=1:jm
            if(double(u(im, j, k) > 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(imm1,j));
                    tf(im,j,k)= t(im,j,k)-2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))*(t(im,j,k)-t(imm1,j,k))-wm*(t(imm1,j,k-1)-t(imm1,j,k+1));
                end
            end
            
            if(double(u(2, j, k) < 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(2,j));
                    tf(1,j,k)=t(1,j,k)-2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))*(t(2,j,k)-t(1,j,k))-wm*(t(2,j,k-1)-t(2,j,k+1));
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
                    tf(i,jm,k)=t(i,jm,k)-2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))*(t(i,jm,k)-t(i,jmm1,k))-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1));                   
                end
            end
            if(double(v(i, 2, k) < 0.e0))
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(i,2));
                    tf(i,1,k)=t(i,1,k)-2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))*(t(i,2,k)-t(i,1,k))-wm*(t(i,2,k-1)-t(i,2,k+1));
                end
            end
        end
    end

	tf = tf .* fsm_3d;
    
	return

