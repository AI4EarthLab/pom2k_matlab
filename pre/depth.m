function depth()
global z zz dz dzz z_3d zz_3d dz_3d dzz_3d im jm kb kl1 kl2;
kdz=[1,1,2,4,8,16,32,64,128,256,512,1024];

for k=2:kl1
    z(k)=z(k-1)+kdz(k-1);
end

delz=z(kl1)-z(kl1-1);

for k=kl1+1:kl2
    z(k)=z(k-1)+delz;
end

for k=kl2+1:kb
    dz(k)=kdz(kb-k+1)*delz/kdz(kb-kl2);
    z(k)=z(k-1)+dz(k);
end

for k=1:kb
    z(k)=-z(k)/z(kb);
end

for k=1:kb-1
    zz(k)=0.5e0*(z(k)+z(k+1));
end

zz(kb)=2.e0*zz(kb-1)-zz(kb-2);

for k=1:kb-1
    dz(k)=z(k)-z(k+1);
    dzz(k)=zz(k)-zz(k+1);
end

dz(kb)=0.e0;    dzz(kb)=0.e0;

for k=1:kb
    z_3d(:,:,k)=repmat(z(k),im,jm);     zz_3d(:,:,k)=repmat(zz(k),im,jm);
    dz_3d(:,:,k)=repmat(dz(k),im,jm);   dzz_3d(:,:,k)=repmat(dzz(k),im,jm);
end
return