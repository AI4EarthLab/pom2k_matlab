function [tf] = advt2(tb,t,tclim,tf,nitera,sw,etb,etf,w,dt,dt_3d,aam,tprni,h,dum_3d,dvm_3d,u,v)
global im kb kbm1 dti2 dx_3d dy_3d;
              
    xmassflux=AXB(dt_3d).*u;        ymassflux=AYB(dt_3d).*v;        eta=etb;
	tb(:,:,kb)=tb(:,:,kbm1);        zwflux=w;                       fbmem=tb;
    tmp1=tf(1,:,:);                 tmpim=tf(im,:,:);               h_3d=repmat(h,1,1,kb);
%
%     Start Smolarkiewicz scheme:
%
for itera=1:nitera
    %     Upwind advection scheme:
    xflux=(xmassflux>0.e0).* xmassflux .* shift(fbmem,1,1) + (xmassflux<0.e0).* xmassflux .* fbmem ;
    yflux=(ymassflux>0.e0).* ymassflux .* shift(fbmem,1,2) + (ymassflux<0.e0).* ymassflux .* fbmem ;   
    zflux=((zwflux < 0.e0).* zwflux .* shift(fbmem,1,3) + (zwflux > 0.e0).* zwflux .* fbmem );
    zflux(:,:,1)=0.e0;
    if(itera==1)
       zflux(:,:,1)=w(:,:,1).*t(:,:,1);
    end
    zflux(:,:,kb)=0.e0; 
    %
    %     Add net advective fluxes and step forward in time:
    %
    tf=(fbmem .* repmat(h+eta,1,1,kb) - dti2.*( DXF(xflux) + DYF(yflux) - DZF(zflux) ) ) ./ repmat(h+etf,1,1,kb) ;
    tf(1,:,:)=tmp1; tf(im,:,:)=tmpim;
    %
    %     calculate antidiffusion velocity:
    %
%    [xmassflux,ymassflux,zwflux,tf] = smol_adif(xmassflux,ymassflux,zwflux,tf,sw,dt);

    [xmassflux,ymassflux,zwflux,tf] = smol_adif(xmassflux.*AXB(dy_3d),ymassflux.*AYB(dx_3d),zwflux,tf,sw,dt);
    xmassflux=xmassflux./AXB(dy_3d);    ymassflux=ymassflux./AYB(dx_3d);    
    eta=etf;       fbmem=tf;
    %     End of Smolarkiewicz scheme
end
%
%     Add horizontal diffusive fluxes:
%
tb=tb-tclim;    xmassflux=AXB(aam);     ymassflux=AYB(aam);

xflux=-xmassflux.*AXB(h_3d).*tprni.*DXB(tb).*dum_3d;
yflux=-ymassflux.*AYB(h_3d).*tprni.*DYB(tb).*dvm_3d;

tf=tf-dti2.* ( DXF(xflux)+ DYF(yflux) ) ./ repmat(h+etf,1,1,kb);
tf(:,:,kb)=0.e0;    tf(im,:,:)=tmpim;

end
