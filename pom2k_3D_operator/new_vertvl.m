function [xflux,yflux,w]=new_vertvl (xflux,yflux,w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  calculates vertical velocity.                       *
% *                                                                    *
% **********************************************************************
load('grid.mat'); load('operator.mat');
dt_3d=repmat(dt,1,1,kb);
w = zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);
%     Reestablish boundary conditions:
% % %       for k=1:kbm1
% % %         for j=2:jm
% % %           for i=2:im
% % %             xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)) ...
% % %                     *(dt(i,j)+dt(i-1,j))*u(i,j,k);
% % %           end
% % %         end
% % %       end
% % %     
% % % 
% % % %
% % %       for k=1:kbm1
% % %         for j=2:jm
% % %           for i=2:im
% % %             yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1))     ... 
% % %                    *(dt(i,j)+dt(i,j-1))*v(i,j,k);
% % %           end
% % %         end
% % %       end
xflux = AXB(dy_3d) .* AXB(dt_3d) .* u; 
yflux = AYB(dx_3d) .* AYB(dt_3d) .* v;

%     NOTE that, if one wishes to include freshwater flux, the
%     surface velocity should be set to vflux(i,j). See also
%     change made to 2-D volume conservation equation which
%     calculates elf.
%
        for j=2:jmm1
          for i=2:imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j));
          end
        end

w(:,:,1)=0.5*(vfluxb+vfluxf);

w0=w;     
%       for k=1:kbm1
%         for j=2:jmm1
%           for i=2:imm1
%             w(i,j,k+1)=w(i,j,k)     ...
%                 +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)     ...  
%                       +yflux(i,j+1,k)-yflux(i,j,k))     ... 
%                        /(dx(i,j)*dy(i,j))     ... 
%                        +(etf(i,j)-etb(i,j))/dti2);
%           end
%         end
%       end

tps=OP_L_XZ*(etf-etb)/dti2;
tps=repmat(tps,1,1,kb);
w=SUM1(w+dz_3d .* ( ( DXF(xflux)+DYF(yflux) )./(dx_3d.*dy_3d)+ tps));

end
