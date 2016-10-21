function getpageprofvel(ij,file,itime)
%
% getpageprofvel:  gets and plots a page of velocity profiles
%
% Usage: getpageprofvel(ij,file,itime)
%
% where: ij ...... matrix of i and j cell indices for profiles
%        file .... the name of the netCDF file
%        itime ... the time index (first output = 1, etc.)
%
% Initial version, JRH 11/12/2001
% itime added to title, JRH 02/01/2002
%
% Number of subplots across page:
%
nxplot=3;
%
% Number of subplots down page:
%
nyplot=4;
%
% Number of subplots:
%
nplot=size(ij,1);
%
% Width of each subplot (normalised):
%
subwidth=0.22;
%
% Height of each subplot (normalised):
%
subheight=0.16;
%
% Horizontal offset of each subplot (normalised):
%
xoff=0.33;
%
% Vertical offset of each subplot (normalised):
%
yoff=0.24;
%
% Left margin (normalised):
%
xmarg=0.1;
%
% Top margin (normalised):
%
ymarg=0.05;
%
% Page height (normalised):
%
pageheight=0.98;
%
set( gcf , 'Position' , [360 187 560 747] );
set( gcf , 'PaperPosition' , [0.5 0.5 7.25 10.5] );
set(gca,'FontSize',16);
%
for iplot = 1:nplot,
%
  iyplot=floor((iplot-1)/nxplot)+1;
  ixplot=iplot-(iyplot-1)*nxplot;
  subplot('position',[(ixplot-1)*xoff+xmarg pageheight-iyplot*yoff+ymarg subwidth subheight ]);
  set(gca,'FontSize',12);
  pbaspect([1 1 1]);
  i = ij(iplot,1);
  j = ij(iplot,2);
  [vel,height,heightlim] = profvel(file,itime,i,j);
  pprofvel(vel,height,heightlim,i,j,itime,iplot,0);
%
end
%
% Final plot for legend only:
%
iplot=nplot+1;
iyplot=floor((iplot-1)/nxplot)+1;
ixplot=iplot-(iyplot-1)*nxplot;
subplot('position',[(ixplot-1)*xoff+xmarg pageheight-iyplot*yoff+ymarg subwidth subheight ]);
set(gca,'FontSize',12);
hold on;
tag=plot(0, 0,'rx-');
set(tag,'Visible','off');
tag=plot(0,0,'b+-');
set(tag,'Visible','off');
legend('Eastward','Northward');
hold off;
set(gca,'Visible','off');
