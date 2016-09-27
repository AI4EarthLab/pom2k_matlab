function [advx,advy]=new_advct(u,v,dx,dy,dt,aam,ub,vb,aru,arv)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates the horizontal portions of momentum      *
% *                advection well in advance of their use in advu and  *
% *                advv so that their vertical integrals (created in   *
% *                the main program) may be used in the external (2-D) *
% *                mode calculation.                                   *
% *                                                                    *
% **********************************************************************
%
%
load('grid.mat');

curv=zeros(im,jm,kb);
advx=zeros(im,jm,kb);
advy=zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);

xflux1=zeros(im,jm,kb);
yflux1=zeros(im,jm,kb);
advx1=zeros(im,jm,kb);
%
for k=1:kb
    curv(:,:,k)= (AYF2_XY(v(:,:,k)) .* DXB2_XY(AXF2_XY(dy)) - AXF2_XY(u(:,:,k)) .* DYB2_XY(AYF2_XY(dx)) ) ./ (dx .* dy);
end

%
%     Calculate x-component of velocity advection:
%
%     Calculate horizontal advective fluxes:
%

for k=1:kbm1
    for j=1:jm
        for i=2:imm1
            xflux(i,j,k)=.125e0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)     ...
                +(dt(i,j)+dt(i-1,j))*u(i,j,k))     ...
                *(u(i+1,j,k)+u(i,j,k));
        end
    end
end

for k=1:kbm1
    for j=2:jm
        for i=2:im
            yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k)     ...
                +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))     ...
                *(u(i,j,k)+u(i,j-1,k));
        end
    end
end

%
%    Add horizontal diffusive fluxes:
%
for k=1:kbm1
    for j=2:jm
        for i=2:imm1
            xflux(i,j,k)=xflux(i,j,k)     ...
                -dt(i,j)*aam(i,j,k)*2.e0     ...
                *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j);
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))     ...
                *(aam(i,j,k)+aam(i-1,j,k)     ...
                +aam(i,j-1,k)+aam(i-1,j-1,k));
            yflux(i,j,k)=yflux(i,j,k)     ...
                -dtaam*((ub(i,j,k)-ub(i,j-1,k))     ...
                /(dy(i,j)+dy(i-1,j)     ...
                +dy(i,j-1)+dy(i-1,j-1)) ...
                +(vb(i,j,k)-vb(i-1,j,k))...
                /(dx(i,j)+dx(i-1,j)     ...
                +dx(i,j-1)+dx(i-1,j-1)));
            %
            xflux(i,j,k)=dy(i,j)*xflux(i,j,k);
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)  ...
               +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k);
        end
    end
end


for k=1:kb
    xflux1(:,:,k) = dy .* (AXF1_XY( AXB1_XY(dt) .* u(:,:,k) ) .* AXF2_XY(u(:,:,k)) ...
                    - dt.*aam(:,:,k)*2.0.*DXF2_XY(ub(:,:,k))./dx );
    tmp1= ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k)) );
    tmp2=AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
                     .*( DYB1_XY(ub(:,:,k)) ./ AYB1_XY(AXB2_XY(dy)) + DXB1_XY(vb(:,:,k)) ./ AYB1_XY(AXB2_XY(dx))) ;
    tmp2(isnan(tmp2))=0;
    tmp2(isinf(tmp2))=0;             
   % yflux1(:,:,k) = AYB2_XY(AXB1_XY(dx)) .* (tmp1-tmp2) ;  
   yflux1(:,:,k) = AYB1_XY(AXB2_XY(dx)) .* (tmp1-tmp2) ;
end


%
%     for horizontal advection:
%

for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)     ...
                +yflux(i,j+1,k)-yflux(i,j,k);
        end
    end
end

for k=1:kbm1
    for j=2:jmm1
        for i=3:imm1
            advx(i,j,k)=advx(i,j,k)  ...
                -aru(i,j)*.25e0     ...
                *(curv(i,j,k)*dt(i,j)     ...
                *(v(i,j+1,k)+v(i,j,k))     ...
                +curv(i-1,j,k)*dt(i-1,j)     ...
                *(v(i-1,j+1,k)+v(i-1,j,k)));
        end
    end
end

for k=1:kb
    xflux1(:,:,k) = dy .* (AXF2_XY( AXB1_XY(dt) .* u(:,:,k) ) .* AXF2_XY(u(:,:,k)) ...
                    - dt.*aam(:,:,k)*2.0.*DXF2_XY(ub(:,:,k))./dx );
    yflux1(:,:,k) = AYB1_XY(AXB1_XY(dx)) .* ( ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k))) ...
                    - AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
                    .* ( DYB1_XY(ub(:,:,k)) ./ AYB1_XY(AXB2_XY(dy)) + DXB1_XY(vb(:,:,k)) ./ AYB1_XY(AXB2_XY(dx)) ) );                
    %set NaN caused by D/A operation to 0.
    yflux1(isnan(yflux1))=0; 
    yflux1(isinf(yflux1))=0;
    %     for horizontal advection:   
    advx1(:,:,k)=DXB1_XY( xflux1(:,:,k) ) + DYF1_XY( yflux1(:,:,k) ) ...
                    - aru .* AXB1_XY(  curv(:,:,k) .* dt .* AYF1_XY(v(:,:,k)) );    
end

%
%-----------------------------------------------------------------------
%
advy=zeros(im,jm,kb);
xflux=zeros(im,jm,kb);
yflux=zeros(im,jm,kb);
%
%     Calculate y-component of velocity advection:
%
%     Calculate horizontal advective fluxes:
%
for k=1:kbm1
    for j=2:jm
        for i=2:im
            xflux(i,j,k)=.125e0*((dt(i,j)+dt(i-1,j))*u(i,j,k)     ...
                +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))     ...
                *(v(i,j,k)+v(i-1,j,k));
        end
    end
end

for k=1:kbm1
    for j=2:jmm1
        for i=1:im
            yflux(i,j,k)=.125e0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)     ...
                +(dt(i,j)+dt(i,j-1))*v(i,j,k))     ...
                *(v(i,j+1,k)+v(i,j,k));
        end
    end
end

%
%    Add horizontal diffusive fluxes:
%
for k=1:kbm1
    for j=2:jmm1
        for i=2:im
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))     ...
                *(aam(i,j,k)+aam(i-1,j,k)     ...
                +aam(i,j-1,k)+aam(i-1,j-1,k));
            xflux(i,j,k)=xflux(i,j,k)     ...
                -dtaam*((ub(i,j,k)-ub(i,j-1,k))     ...
                /(dy(i,j)+dy(i-1,j)     ...
                +dy(i,j-1)+dy(i-1,j-1))     ...
                +(vb(i,j,k)-vb(i-1,j,k))     ...
                /(dx(i,j)+dx(i-1,j)     ...
                +dx(i,j-1)+dx(i-1,j-1)));
            yflux(i,j,k)=yflux(i,j,k)     ...
                -dt(i,j)*aam(i,j,k)*2.e0     ...
                *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j);
            %
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)     ...
                +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k);
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k);
        end
    end
end

%
%     for horizontal advection:
%

for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)     ...
                +yflux(i,j,k)-yflux(i,j-1,k);
        end
    end
end

for k=1:kbm1
    for j=3:jmm1
        for i=2:imm1
            advy(i,j,k)=advy(i,j,k)     ...
                +arv(i,j)*.25e0     ...
                *(curv(i,j,k)*dt(i,j)     ...
                *(u(i+1,j,k)+u(i,j,k))     ...
                +curv(i,j-1,k)*dt(i,j-1)     ...
                *(u(i+1,j-1,k)+u(i,j-1,k)));
        end
    end
end
