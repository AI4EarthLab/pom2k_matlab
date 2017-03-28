    function [tf,t,tb]=cal_t(tb,t,dt,dt_3d,u,v,aam,w,etb,etf,wtsurf,tsurf,nbct,swrad,kh)
    global kb dum_3d dvm_3d dti2 tprni h_3d im jm dz_3d dzz_3d kbm1 umol kbm2 dz gs h tbe tbw tbn tbs;
    %Explicit solution
    tf=((h_3d+repmat(etb,1,1,kb)).*tb - dti2 .* (DXF( AXB(dt_3d).*AXB(t).*u-AXB(aam).*AXB(h_3d)*tprni.*DXB(tb).*dum_3d ) ...
        + DYF( AYB(dt_3d).*AYB(t).*v-AYB(aam).*AYB(h_3d)*tprni.*DYB(tb).*dvm_3d )-DZF( AZB(t).*w ))) ./((h_3d+repmat(etf,1,1,kb))) ;
    %-----------------------------------------------------------------
    %Implicit solution
    r=[0.58,0.62,0.67,0.77,0.78];
    ad1=[0.35,0.60,1.0,1.5,1.4];
    ad2=[0.23,20.0,17.0,14.0,7.90];

    rad=create_field(zeros(im,jm,kb),gs,7);
    a = create_field(zeros(im,jm,kb),gs,7);
    c = create_field(zeros(im,jm,kb),gs,7);
    dh = h+etf;     ee=zeros(im,jm,kb);     gg=zeros(im,jm,kb);
    dh_3d=repmat(dh,1,1,kb);

    a(:,:,1:kbm2)=-dti2*(kh(:,:,2:kbm1)+umol);
    a=DIVISION(a,dz_3d.*dzz_3d.*dh_3d.*dh_3d);
    c(:,:,2:kbm1)=dzz_3d(:,:,1:kbm2);
    c=DIVISION(-dti2*(kh+umol),dz_3d.*c.*dh_3d.*dh_3d);
    if(nbct==2||nbct==4)
        rad=repmat(swrad,1,1,kb).*(r(ntp)*exp(z_3d.*dh_3d/ad1(ntp))+(1.e0-r(ntp))*exp(z_3d.*dh_3d/ad2(ntp)));
        rad(:,:,kb)=0.0;
    end
    if(nbct==1)
        ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
        gg(:,:,1)= (-dti2*wtsurf./(-dz(1)*dh)-tf(:,:,1))./(a(:,:,1)-1.e0);
    elseif(nbct==2)
        ee(:,:,1)=a(:,:,1)./(a(:,:,1)-1.e0);
        gg(:,:,1)= (dti2*(wtsurf+rad(:,:,1)-rad(:,:,2))./(dz(1)*dh)-tf(:,:,1))./(a(:,:,1)-1.e0);
    elseif(nbct==3 || nbct==4)
        ee(:,:,1)=0.e0;
        gg(:,:,1)=tsurf;
    end
    for k=2:kbm2
        gg(:,:,k)=1.e0./(a(:,:,k)+c(:,:,k).*(1.e0-ee(:,:,k-1))-1.e0);
        ee(:,:,k)=a(:,:,k).*gg(:,:,k);
        gg(:,:,k)=(c(:,:,k).*gg(:,:,k-1)-tf(:,:,k)+dti2*(rad(:,:,k)-rad(:,:,k+1))./(dh*dz(k))).*gg(:,:,k);
    end
    tf(:,:,kbm1)=(c(:,:,kbm1).*gg(:,:,kbm2)-tf(:,:,kbm1)+dti2*(rad(:,:,kbm1)-rad(:,:,kb)) ...
        ./(dh*dz(kbm1)))./(c(:,:,kbm1).*(1.e0-ee(:,:,kbm2))-1.e0);
    for k=kbm2:-1:1
        tf(:,:,k)=ee(:,:,k).*tf(:,:,k+1)+gg(:,:,k);
    end

%     [tf] = bcond4_t(tf,u,v,w,t,dt,tbe,tbw,tbn,tbs);
    [tf] = bcond4(tf,u,v,w,t,dt,tbe,tbw,tbn,tbs);
    [t,~,tb,~]=smoth_update(tf,0,t,0,tb,0);
    return