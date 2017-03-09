function [e_atmos,vfluxf,wtsurf,swrad,wssurf,advx,advy,drhox,drhoy,aam,adx2d,ady2d,drx2d,dry2d,aam2d,advua,advva,w,egf,utf,vtf]= ...
         set_lateralboundary(e_atmos,vfluxf,wtsurf,swrad,wssurf,u,ub,v,vb,dt_3d,aam,rho,rmean,ramp,uab,vab,ua,va,d,t,s,w,el)
global im jm kb tbias sbias dz_3d dx_3d dy_3d horcon aam_init ispi isp2i;
e_atmos(:,:)=0.e0;  vfluxf(:,:)=0.e0;    wtsurf(:,:)=0.0;     
swrad(:,:)=0.0;     satm=0.e0;
tatm=zeros(im,jm);
w(:,:,1)=vfluxf;
tatm(:,:)=t(:,:,1)+tbias ;
wtsurf=wtsurf+vfluxf.*(tatm-t(:,:,1)-tbias);
wssurf=  vfluxf.*(satm-s(:,:,1)-sbias)  ;
    
        [advx,  advy]   = advct (u,v,dt_3d,aam,ub,vb);        
        [drhox,drhoy]   = baropg(rho, rmean, dt_3d, ramp);  
        
        aam=horcon .* dx_3d .* dy_3d .*sqrt( DXF(u).^2 + DYF(v).^2 +0.5*( DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))) ).^2 );     
        aam(:,1,:)=aam_init;
        aam(:,jm,:)=aam_init;
        aam(1,:,:)=aam_init;
        aam(im,:,:)=aam_init;
        aam(:,:,kb)=aam_init;
                
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
return
end