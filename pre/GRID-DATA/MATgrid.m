%
% MATLAB code to plot output from GRID.f  (T.Ezer, 3/4/2003)
%
% ------- READ GRID PARAMETERS -------
 load ijk.dat
 im=ijk(1);, jm=ijk(2);, kb=ijk(3);
 for k=1:kb, zz(k)=ijk(3+k);, end
 t(im*jm,kb)=0.;
%
% ------- READ BOUNDARY VALUES SET in gridborder -------
 load bnd.dat
 xb=bnd(:,1);, yb=bnd(:,2);
%
% ------- READ GRID, BATHYMETRY & INITIAL TEMP. ------
%
 load plt.dat
 lon=plt(:,3);, lat=plt(:,4);,h=plt(:,5);, for k=1:kb, t(:,k)=plt(:,5+k);, end
%
% ------- READ DX and DY (in km)                ------
%
 load dxy.dat
 dxx=dxy(:,5);, dyy=dxy(:,6);
%
% ------- READ WIND VELOCITY ------
%
 load wnd.dat
 uu=wnd(:,5);, vv=wnd(:,6);
%
  n=0;
 for i=1:im, for j=1:jm, n=n+1;
  x(i,j)=lon(n);, y(i,j)=lat(n);,d(i,j)=h(n);, d1(i,j)=-1*h(n);
  u(i,j)=uu(n);, v(i,j)=vv(n);, dx(i,j)=dxx(n);, dy(i,j)=dyy(n); 
  ts(i,j)=t(n,1);, tb(i,j)=t(n,kb);
 if d1(i,j) == -1, d1(i,j)=max(h)/64;, end   %--- for land mask ---
% -- extract section
 if j == 20
  for k=1:kb, tt(i,k)=t(n,k);, end
 end
 end, end

% ------- PLOT GRID on top of BATHYMETRY DATA used as input ------
figure
% ---- Small domain (-79W,-70W) & (31N,40N) ---- 
 load TOPO.dat, e=TOPO;, xy=[-79 -70 31 40];, imt=109;, jmt=109;
 for i=1:imt, xt(i)=xy(1)+(i-1)/12;, end
 for j=1:jmt, yt(j)=xy(3)+(j-1)/12;, end

 e=flipud(e);,  %e=min(e,3000);, e=max(e,-6000);
 pcolor(xt,yt,e), axis equal, shading flat, colorbar, axis(xy), hold on
 contour(xt,yt,e,[0 0],'w')     % contour of coastline
 caxis([-5000 1000])
 title('ETOPO5 Elevation Data')
 xlabel('Longitude')
 ylabel('Latitude')
% now add model grid on top of topography
  for i=1:im, plot(x(i,:),y(i,:),'k'), end
  for j=1:jm, plot(x(:,j),y(:,j),'k'), end

% ------- PLOT GRID & MODEL BATHYMETRY ------
figure
%pcolor(x,y,d), colorbar, grid on, hold on
% --- to mask land in white, modify colormap
 c=colormap;, c(64,:)=1;
 pcolor(x,y,d1), colorbar, colormap(c), grid on, hold on
 title('Model Grid and Bottom Topography')
 xlabel('Longitude')
 ylabel('Latitude')
% show boundary points from gridborder
  plot(xb,yb,'d','MarkerSize',[12],'MarkerFaceColor',['g'])
% --- save postscript
% print -dpsc PLOTgrid.ps

% ------- PLOT GRID SIZE DX & DY             ---
figure
 pcolor(x,y,dx), shading('flat'), colorbar, grid on, hold on
 [cs,h]=contour(x,y,dx,'k-');, clabel(cs,h)
 title(' DX (km) ')
 xlabel('Longitude')
 ylabel('Latitude')
figure
 pcolor(x,y,dy), shading('flat'), colorbar, grid on, hold on
 [cs,h]=contour(x,y,dy,'k-');, clabel(cs,h)
 title(' DY (km) ')
 xlabel('Longitude')
 ylabel('Latitude')

% ------- PLOT SURFACE TEMPERATURE AND WIND  ---
figure
 pcolor(x,y,ts), shading('interp'), colorbar, grid on
 hold on, quiver(x,y,u,v)
 title('Surface Temperature and Wind ')
 xlabel('Longitude')
 ylabel('Latitude')
% --- save postscript
% print -dpsc PLOTsurf.ps

% ------- PLOT TEMPERATURE ON BOTTOM SIGMA LEVEL ---
figure
 pcolor(x,y,tb), shading('interp'), colorbar, grid on
 title('Bottom Temperature ')
 xlabel('Longitude')
 ylabel('Latitude')
% --- save postscript
% print -dpsc PLOTbot.ps

% ------- PLOT TEMPERATURE CROSS-SECTION ---
figure
 for i=1:im, for k=1:kb
  az(i,k)=d(i,20)*zz(k);
  xx(i,k)=x(i,20);
 end, end
 pcolor(xx,az,tt), colorbar
 title('Temperature Cross-Section (J=20)')
 xlabel('Longitude')
 ylabel('Depth (m)')
% --- save postscript
% print -dpsc PLOTsec.ps

