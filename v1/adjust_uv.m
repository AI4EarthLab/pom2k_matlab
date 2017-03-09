function [u,v] = adjust_uv(u,utb,utf,v,vtb,vtf,dt_3d)
global kb dz_3d
            u = u - repmat(sum(u.*dz_3d, 3),1,1,kb) + DIVISION(repmat(utb+utf,1,1,kb),2.0 * AXB(dt_3d));            
            v = v - repmat(sum(v.*dz_3d, 3),1,1,kb) + DIVISION(repmat(vtb+vtf,1,1,kb),2.0 * AYB(dt_3d));
            u(:,:,kb)=0.0;    v(:,:,kb)=0.0;
return