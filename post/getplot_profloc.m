function getplot_profloc(ij,file)
%
% getplot_profloc: gets and plots locations of cells for velocity profiles
%
% Usage: getplot_profloc(ij,file)
%
% where: ij ...... matrix of ix and iy cell indices for profiles
%        file .... name of input NetCDF file
%
% Initial version, JRH 11/12/2001
%
hold on;
%
if (nargin ~= 2)
  help getplot_profloc;
  return
end
%
% Turn off warnings from netCDF:
%
ncmex('setopts',0);
%
% Open netCDF file:
%
ncid=ncmex('open',file,'nowrite');
%
if(ncid==-1)
  disp(['File ' file ' not found'])
  return
end
%
[name,nx]=ncmex('diminq',ncid,'x');
[name,ny]=ncmex('diminq',ncid,'y');
%
% Get horizontal coordinates:
%
east_e=ncmex('varget',ncid,'east_e',[0 0],[-1 -1],1);
north_e=ncmex('varget',ncid,'north_e',[0 0],[-1 -1],1);
east_c=ncmex('varget',ncid,'east_c',[0 0],[-1 -1],1);
north_c=ncmex('varget',ncid,'north_c',[0 0],[-1 -1],1);
%
% Extrapolate end values:
%
right=2*east_c(nx,:)-east_c(nx-1,:);
top=2*east_c(:,ny)-east_c(:,ny-1);
top_right=2*east_c(nx,ny)-east_c(nx-1,ny-1);
east_c=[ east_c top ; right top_right ];
%
right=2*north_c(nx,:)-north_c(nx-1,:);
top=2*north_c(:,ny)-north_c(:,ny-1);
top_right=2*north_c(nx,ny)-north_c(nx-1,ny-1);
north_c=[ north_c top ; right top_right ];
%
% Get free surface mask:
%
fsm=ncmex('varget',ncid,'fsm',[0 0],[-1 -1],1); 
fsm=[ fsm zeros(nx,1); zeros(1,ny+1) ];
%
% Color in mask:
%
pcolor(east_c,north_c,-fsm);colormap('jet');shading('flat');caxis([-2 .8]);
%
% Number of subplots:
nplot=size(ij,1);
%
for iplot = 1:nplot,  
  xp=east_e(ij(iplot,1),ij(iplot,2));
  yp=north_e(ij(iplot,1),ij(iplot,2));
  plot(xp,yp,'r+ ');
  str=sprintf(' %s',char(iplot+64));
  text(xp,yp,str);
end
%
title('Locations of Current Profiles');
xlabel('metres');
ylabel('metres');
set ( gca, 'DataAspectRatio', [1 1 1] );
axis image;
box on;
hold off;
%
% Close file:
%
ncmex('close',ncid);
