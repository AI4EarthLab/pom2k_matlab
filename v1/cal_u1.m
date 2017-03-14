    function [uf,wubot] =cal_u1(advx,dt_3d,e_atmos,drhox,ub,u,v,w,egf,egb,etf,etb,km,vb,wusurf,cbc)
    global im kb dti2 grav h_3d jm dz_3d dzz_3d kbm1 umol kbm2 dz dum_3d gs h cor ee gg;
    % Explicit solution
    uf=DIVISION( (AXB(repmat(etb,1,1,kb)+h_3d).*ub -dti2*( advx + drhox - AXB( repmat(cor,1,1,kb) .*dt_3d.*AYF(v) )...
        +grav*AXB(dt_3d).*( DXB( repmat(egf+egb,1,1,kb) )+DXB( repmat(e_atmos,1,1,kb) )*2.0 )/2.e0-DZF(AXB(w) .* AZB(u)))), ...
        AXB( repmat(etf,1,1,kb)+h_3d ) );
    uf(:,:,kb)=0.e0;  bond=AXB(w) .* AZB(u); uf(im,:,:) = bond(im,:,:) ;   %add by hx
    %Implicit solution
    dh = AXB(h+etf);dh(1,:)=1.e0;dh(:,1)=1.e0;
    dh_3d=repmat(dh,1,1,kb);
    c=AXB(km);
    a = create_field(zeros(im,jm,kb),gs,6);
    a(:,:,1:kbm2)=-dti2*(c(:,:,2:kbm1)+umol);
    d = create_field(zeros(im,jm,kb),gs,2);
    d(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
    c=DIVISION(-dti2*(c+umol),dz_3d.*d.*dh_3d.*dh_3d);
    ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
    gg(:,:,1)=(-dti2*wusurf./(-dz(1)*dh)-uf(:,:,1))./(a(:,:,1)-1.e0);
    for k=2:kbm2
        gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
        ee(:,:,k)=a(:,:,k).*gg(:,:,k);
        gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-uf(:,:,k)).*gg(:,:,k);
    end
    tps=AXB(cbc).* sqrt( ub(:,:,kbm1).^2 + AXB( AYF( vb(:,:,kbm1) ) ).^2 );
    uf(:,:,kbm1)=DIVISION(c(:,:,kbm1).*gg(:,:,kbm2)-uf(:,:,kbm1), ...
        tps*dti2./(-dz(kbm1)*dh)-1.e0-(ee(:,:,kbm2)-1.e0).*c(:,:,kbm1));
    for k=kbm2:-1:1
        uf(:,:,k)=(ee(:,:,k).*uf(:,:,k+1)+gg(:,:,k));
    end
    uf=uf.*dum_3d;
    wubot=-tps.*uf(:,:,kbm1);

    [uf] = bcond3_u(uf,u);

    return

