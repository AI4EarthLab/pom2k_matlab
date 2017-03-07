            tps=sum(u(:,:,1:kbm1).*dz_3d(:,:,1:kbm1), 3);
            tps1=sum(u.*dz_3d, 3);            
            
            utb_3d = repmat(utb, 1, 1, kbm1);
            utf_3d = repmat(utf, 1, 1, kbm1);
            tps_3d = repmat(tps, 1, 1, kbm1);
            dt_3d1 = repmat(dt, 1, 1, kbm1);
            dt_axb = 2.0 * AXB(dt_3d1);
%             u(:,:,1:kbm1) = u(:,:,1:kbm1)-tps_3d + (utb_3d+utf_3d) ./ dt_axb;
            u(:,:,1:kbm1) = u(:,:,1:kbm1)-tps_3d + DIVISION(utb_3d+utf_3d, dt_axb);            
            tps = sum(v(:,:,1:kbm1) .* dz_3d(:,:,1:kbm1), 3);
            
           
            vtb_3d = repmat(vtb, 1, 1, kbm1);
            vtf_3d = repmat(vtf, 1, 1, kbm1);
            tps_3d = repmat(tps, 1, 1, kbm1);
            dt_ayb = 2.0 * AYB(dt_3d1);
            v(:,:,1:kbm1) = v(:,:,1:kbm1) - tps_3d + (vtb_3d+vtf_3d) ./ dt_ayb;