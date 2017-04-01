    function [vf,wvbot] = internal_v(advy,dt_3d,e_atmos,drhoy,vb,u,v,w,egf,egb,etf,etb,km,ub,wvsurf,cbc)
    global kb dti2 grav h_3d im jm dz_3d dzz_3d kbm1 umol kbm2 dz dvm_3d gs h cor ee gg;
    % Explicit solution
    vf = ( AYB( repmat(etb,1,1,kb) +h_3d).*vb-dti2*( advy+drhoy+AYB( repmat(cor,1,1,kb) .*dt_3d.*AXF(u) )...
        +grav*AYB(dt_3d).*( DYB( repmat(egf+egb,1,1,kb) )+DYB( repmat(e_atmos,1,1,kb) )*2.0 )/2.e0-DZF(AYB(w) .* AZB(v)))) ...
        ./AYB( repmat(etf,1,1,kb) +h_3d);
    vf(:,:,kb)=0.e0;
    %Implicit solution
    dh = AYB(h+etf);dh(1,:)=1.e0;dh(:,1)=1.e0;
    dh_3d=repmat(dh,1,1,kb);
    c=AYB(km);
    a = create_field(zeros(im,jm,kb),gs,1);
    a(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
    d = create_field(zeros(im,jm,kb),gs,1);
    d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
    c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);
    ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
    gg(:,:,1)=(-dti2*wvsurf./(-dz(1)*dh)-vf(:,:,1))./(a(:,:,1)-1.e0);
    for k=2:kbm2
        gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
        ee(:,:,k)=a(:,:,k).*gg(:,:,k);
        gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-vf(:,:,k)).*gg(:,:,k);
    end
    tps=AYB(cbc).* sqrt( AYB( AXF( ub(:,:,kbm1) ) ).^2 + vb(:,:,kbm1).^2 );
    vf(:,:,kbm1)=DIVISION(c(:,:,kbm1).*gg(:,:,kbm2)-vf(:,:,kbm1), ...
        tps*dti2./(-dz(kbm1)*dh)-1.e0-(ee(:,:,kbm2)-1.e0).*c(:,:,kbm1));
    for k=kbm2:-1:1
        vf(:,:,k)=(ee(:,:,k).*vf(:,:,k+1)+gg(:,:,k));
    end


    vf=vf.*dvm_3d;
    wvbot=-tps .* vf(:,:,kbm1);
    [vf] = bcond3_v(vf, v);
    return