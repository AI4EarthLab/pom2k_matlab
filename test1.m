load('grid.mat');

yflux=zeros(im,jm,kb);
yflux1=zeros(im,jm,kb);


for k=1:kbm1
    for j=2:jm
        for i=2:im
            yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k)     ...
                +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))     ...
                *(u(i,j,k)+u(i,j-1,k));
        end
    end
end

for k=1:kb
    yflux1(:,:,k) =  ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k)) ) ;
end


tmp=zeros(im,jm,kb);
for k=1:kbm1
    for j=2:jm
        for i=2:imm1
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))     ...
                *(aam(i,j,k)+aam(i-1,j,k)     ...
                +aam(i,j-1,k)+aam(i-1,j-1,k));
            yflux(i,j,k)= yflux(i,j,k)...
                -dtaam*((ub(i,j,k)-ub(i,j-1,k))     ...
                /(dy(i,j)+dy(i-1,j)     ...
                +dy(i,j-1)+dy(i-1,j-1)) ...
                +(vb(i,j,k)-vb(i-1,j,k))...
                /(dx(i,j)+dx(i-1,j)     ...
                +dx(i,j-1)+dx(i-1,j-1)));
            tmp(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)  ...
                +dx(i,j-1)+dx(i-1,j-1));
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)  ...
                +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k);
        end
    end
end

for k=1:kb
% % %     yflux1(:,:,k) =  ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k)) )...
% % %                      - AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
% % %                      .*( DYB1_XY(ub(:,:,k)) ./ AYB1_XY(AXB2_XY(dy)) + DXB1_XY(vb(:,:,k)) ./ AYB1_XY(AXB2_XY(dx))) ;                
    tmp1= ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k)) );
    tmp2=AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
                     .*( DYB1_XY(ub(:,:,k)) ./ AYB1_XY(AXB2_XY(dy)) + DXB1_XY(vb(:,:,k)) ./ AYB1_XY(AXB2_XY(dx))) ;
    tmp2(isnan(tmp2))=0;
    tmp2(isinf(tmp2))=0;             
   % yflux1(:,:,k) = AYB2_XY(AXB1_XY(dx)) .* (tmp1-tmp2) ;  
    tmp3(:,:,k)= AYB1_XY(AXB2_XY(dx)) ;
    yflux1(:,:,k) =tmp3(:,:,k) .* (tmp1-tmp2) ;
end
% % % 
% % % % for k=1:kb
% % % %     yflux1(:,:,k) =  ( AXB1_XY( AYB1_XY(dt) .* v(:,:,k) ) .* AYB1_XY(u(:,:,k)) ) ...
% % % %                      - AYB1_XY(AXB1_XY(dt)) .* AYB1_XY(AXB1_XY(aam(:,:,k))) ...
% % % %                      .*( DYB1_XY(ub(:,:,k)) ./ AYB1_XY(AXB2_XY(dy)) + DXB1_XY(vb(:,:,k)) ./ AYB1_XY(AXB2_XY(dx))) ;                
% % % %     %set NaN caused by D/A operation to 0.
% % % %     yflux1(isnan(yflux1))=0;
% % % %     yflux1(isinf(yflux1))=0;
% % % % end
