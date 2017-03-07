            tps = sum((uf+ub-2.e0*u).*dz_3d, 3);
            
            tps_3d = repmat(tps, 1, 1, kbm1);
            u(:,:,1:kbm1)=u(:,:,1:kbm1)+  .5e0*smoth*(uf(:,:,1:kbm1)+ub(:,:,1:kbm1)-2.e0*u(:,:,1:kbm1)-tps_3d);

            tps = sum((vf+vb-2.e0*v).*dz_3d, 3);
            
            tps_3d = repmat(tps,1,1,kbm1);
            v(:,:,1:kbm1) = v(:,:,1:kbm1)+ .5e0*smoth*(vf(:,:,1:kbm1)+vb(:,:,1:kbm1)-2.e0*v(:,:,1:kbm1)-tps_3d);
            
            ub=u;
            u=uf;
            vb=v;
            v=vf;