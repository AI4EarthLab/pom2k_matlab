function [tb,tf] = advt2(tb,t,tclim,tf,nitera,sw,etb,etf,w,art,dt,dt_3d,aam,tprni,h,dum_3d,dvm_3d,u,v,aru,arv)
global im jm kb kbm1 dti2;

    fbmen=zeros(im,jm,kb);          h_3d=repmat(h,1,1,kb);              
    art_3d=repmat(art,1,1,kb);      dtb_3d=repmat(h+etb,1,1,kb);    dtf_3d=repmat(h+etf,1,1,kb);

    xmassflux=AXB(dt_3d).*u;
    ymassflux=AYB(dt_3d).*v;

	tb(:,:,kb)=tb(:,:,kbm1);            
    zwflux=w;
    fbmem=tb;
%
%     Start Smolarkiewicz scheme:
%
for itera=1:nitera
    %     Upwind advection scheme:
    
    xflux=(xmassflux>0.e0).* xmassflux .* shift(fbmem,1,1) + (xmassflux<0.e0).* xmassflux .* fbmem ;
    yflux=(ymassflux>0.e0).* ymassflux .* shift(fbmem,1,2) + (ymassflux<0.e0).* ymassflux .* fbmem ;


% 	zflux(:,:,kb)=0.e0;

    %
%     for k=2:kbm1
%         for j=2:jmm1
%             for i=2:imm1
%                 zflux(i,j,k)=0.5e0 *((zwflux(i,j,k)+abs(zwflux(i,j,k))) *fbmem(i,j,k)+     ...
%                     (zwflux(i,j,k)-abs(zwflux(i,j,k)))*fbmem(i,j,k-1));
%             end
%         end
%     end
%     zflux=zflux.*art_3d;
    
    zflux=((zwflux < 0.e0).* zwflux .* shift(fbmem,1,3) + (zwflux > 0.e0).* zwflux .* fbmem ) .*art_3d ;
    zflux(:,:,1)=0.e0;
    if(itera==1)
       zflux(:,:,1)=w(:,:,1).*t(:,:,1).*art(:,:);
    end
    zflux(:,:,kb)=0.e0;
    %
    %     Add net advective fluxes and step forward in time:
    %
%     for k=1:kbm1
%         for j=2:jmm1
%             for i=2:imm1
% %                 ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k)+(zflux(i,j,k)-zflux(i,j,k+1))/dz(k);
%                 ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)-dti2* ...
%                     (xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k)+(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)) )/((h(i,j)+etf(i,j))*art(i,j));
%             end
%         end
%     end
    
    tf=(fbmen.*dtb_3d.*art_3d- dti2.*( DXF(xflux) + DYF(yflux) - DZF(zflux) ) ) ./ dtf_3d ;
    %
    %     %alculate antidiffusion velocity:
    %
  
   [xmassflux,ymassflux,zwflux,tf] = smol_adif(xmassflux,ymassflux,zwflux,tf,sw,aru,arv,dt);
    %
    fbmem=tf;
    %     End of Smolarkiewicz scheme
end
%
%     Add horizontal diffusive fluxes:
%
tb=tb-tclim;    xmassflux=AXB(aam);     ymassflux=AYB(aam);
%
% for k=1:kbm1
%     for j=2:jm
%         for i=2:im
%             xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni     ...
%                 *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)     ...
%                 *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j));
%             yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni     ...
%                 *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)     ...
%                 *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1));
%         end
%     end
% end
xflux=-xmassflux.*AXB(h_3d).*tprni.*DXB(tb).*dum_3d;
yflux=-ymassflux.*AYB(h_3d).*tprni.*DYB(tb).*dvm_3d;
tb=tb+tclim;
%
%     Add net horizontal fluxes and step forward in time:
%
% for k=1:kbm1
%     for j=2:jmm1
%         for i=2:imm1
%             ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)     ...
%                 +yflux(i,j+1,k)-yflux(i,j,k)) /((h(i,j)+etf(i,j))*art(i,j));
%         end
%     end
% end

tf=tf-dti2.* ( DXF(xflux)+ DYF(yflux) ) ./ dtf_3d;

end
