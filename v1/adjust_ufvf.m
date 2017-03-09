function [u,v,ub,vb] = adjust_ufvf(u,ub,uf,v,vb,vf)
global kb dz_3d smoth
    u=u+.5e0*smoth*(uf+ub-2.e0*u-repmat(sum((uf+ub-2.e0*u).*dz_3d, 3),1,1,kb)); 
    u(:,:,kb)=0.0;

    v = v+ .5e0*smoth*(vf+vb-2.e0*v-repmat(sum((vf+vb-2.e0*v).*dz_3d, 3),1,1,kb));
    v(:,:,kb)=0.0;
    ub=u;   u=uf;   vb=v;   v=vf;
return