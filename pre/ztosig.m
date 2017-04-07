function [tb,sb]=ztosig(t,s,h,zz,z2)
global x_d y_d;
% Construct interpolation grid
[im,jm,ks]=size(t);
[kb]=size(zz,1);
x_d=linspace(0,im-1,im);
y_d=linspace(0,jm-1,jm);
z_d=linspace(0,kb-1,kb);
[Y,X,Z]=meshgrid(y_d,x_d,z2'); %z
[Y1,X1,Z1]=meshgrid(y_d,x_d,z_d); % sigm
%Construct vertical grids in sigma coordinate mapping with Z.
for i=1:im
    for j=1:jm
        zzh(i,j,:)=-zz.*h(i,j)      ;
    end
end
% Temperature and salt interpolation:
tb=interp3(Y,X,Z,t,Y,X,zzh,'neareast');
sb=interp3(Y,X,Z,s,Y1,X1,zzh,'neareast');
end