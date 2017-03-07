        [advx,  advy]   = advct (u,v,dt_3d,aam,ub,vb);        
        [drhox,drhoy]   = baropg(rho, rmean, dt_3d, ramp);  
        
        set_aam();
                
        adx2d = sum(advx.*dz_3d, 3);
        ady2d = sum(advy.*dz_3d, 3);
        drx2d = sum(drhox.*dz_3d, 3);
        dry2d = sum(drhoy.*dz_3d, 3);
        aam2d = sum(aam.*dz_3d, 3);
        
        [advua,advva]   = advave(aam2d,uab,vab,ua,va,d);
        adx2d=adx2d-advua;
        ady2d=ady2d-advva;