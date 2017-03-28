function [adx2d,ady2d,drx2d,dry2d,aam2d,advua,advva,egf,utf,vtf] = mode_interaction(advx,advy,drhox,drhoy,aam,uab,vab,ua,va,el,d)
% form vertical averages of 3-D fields for use in external (2-D) mode
global dz_3d ispi isp2i;

        adx2d = sum(advx.*dz_3d, 3);
        ady2d = sum(advy.*dz_3d, 3);
        drx2d = sum(drhox.*dz_3d, 3);
        dry2d = sum(drhoy.*dz_3d, 3);
        aam2d = sum(aam.*dz_3d, 3);
        
        [advua,advva]   = advave(aam2d,uab,vab,ua,va,d);
        adx2d=adx2d-advua;
        ady2d=ady2d-advva;
        
        egf=el*ispi;
        utf=ua .* 2.0 .* AXB(d) .* isp2i;
        vtf=va .* 2.0 .* AYB(d) .* isp2i;
end