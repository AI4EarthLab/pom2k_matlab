%
% plots topographic elevation map from the ETOPO5 data set
%
% ---- Large domain (-100W,-50W) & (0N,50N) ---- 
%load ETOPO5.dat, e=ETOPO5;, xy=[-100 -50 0 50];, im=601;, jm=601;

% ---- Small domain (-79W,-70W) & (31N,40N) ---- 
 load TOPO.dat, e=TOPO;, xy=[-79 -70 31 40];, im=109;, jm=109;

 for i=1:im, x(i)=xy(1)+(i-1)/12;, end
 for j=1:jm, y(j)=xy(3)+(j-1)/12;, end

 e=flipud(e);,  e=min(e,3000);, e=max(e,-6000); 
 pcolor(x,y,e), axis equal, shading flat, hold on
 contour(x,y,e,[0 0],'b')
 caxis([-6000 3000]), colorbar
 axis(xy)

 title('ETOPO5 Elevation Data')
 xlabel('Longitude')
 ylabel('Latitude')
