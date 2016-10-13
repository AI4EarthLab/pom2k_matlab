clear all;

load('grid.mat');
load('test.mat');

aam=zeros(im,jm,kb)+500;
aam1=aam;

A1=zeros(im,jm,kb);A2=A1;A3=A1;A4=A1;
B1=A1;B2=A1;B3=A1;B4=A1;

for k=1:kbm1
    for j=2:jmm1
        for i=2:imm1
            A1(i,j,k)=((u(i+1,j,k)-u(i,j,k))/dx(i,j))^2;
            A2(i,j,k)=((v(i,j+1,k)-v(i,j,k))/dy(i,j))^2 ;
            A3(i,j,k)=0.25*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))/dy(i,j);
            A4(i,j,k)=0.25*(v(i+1,j,k)+v(i+1,j+1,k)-v(i-1,j,k)-v(i-1,j+1,k))/dx(i,j);
            aam(i,j,k)=horcon*dx(i,j)*dy(i,j) ...
                *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))^2     ...
                +((v(i,j+1,k)-v(i,j,k))/dy(i,j))^2     ...
                +0.5*(0.25*(u(i,j+1,k)+u(i+1,j+1,k)     ...
                -u(i,j-1,k)-u(i+1,j-1,k))     ...
                /dy(i,j)     ...
                +0.25*(v(i+1,j,k)+v(i+1,j+1,k)     ...
                -v(i-1,j,k)-v(i-1,j+1,k))     ...
                /dx(i,j))^2);
        end
    end
end

save('test.mat','horcon','dx','dy','u','v');

for k=1:kbm1
    B1(:,:,k)=(DXF2_XY(u(:,:,k))./dx(:,:)).^2;
    B2(:,:,k)=(DYF2_XY(v(:,:,k))./dy(:,:)).^2;
    B3(:,:,k)=DYB2_XY(AYF1_XY(AXF1_XY(u(:,:,k))))./dy(:,:);
    B4(:,:,k)=DXB2_XY(AXF1_XY(AYF1_XY(v(:,:,k))))./dx(:,:);
    
    aam1(:,:,k)=horcon.*dx(:,:).*dy(:,:)...
                .*sqrt( (DXF2_XY(u(:,:,k))./dx).^2 + (DYF2_XY(v(:,:,k))./dy).^2    ...
                +0.5*( DYB2_XY(AYF1_XY(AXF1_XY(u(:,:,k))))./dy + DXB2_XY(AXF1_XY(AYF1_XY(v(:,:,k))))./dx).^2 );
end  
