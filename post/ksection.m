function [var,east_c,north_c]=ksection(file,vname,itime,kindex)
%
% ksection:        returns a "horizontal" section of a 2-D or 3-D variable 
%                  along k=kindex from a *.nc file from pom2k
%
% Usage: [var,east_c,north_c]=ksection(file,vname,itime,kindex)
%
% where: var ..... the 2-D section
%        east_c .. horizontal coordinates of corners of model "cell" (metres)
%        north_c . horizontal coordinates of corners of model "cell" (metres)
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%                  (use negative value if not required)
%        kindex .. k-index across which the section is taken (increasing 
%                  downwards from 1 at the surface)
%                  (use negative value if not required)
%    
% Examples: 
%
%        [s,east_c,north_c]=ksection('file.nc','s',2,3)
%           gets the salinity for time 2 at k=3 from file file.nc
%
%        [elb,east_c,north_c]=ksection('file.nc','elb',2,-1)
%           gets the suface elevation for time 2 from file file.nc
%
%        [h,east_c,north_c]=ksection('file.nc','h',-1,-1)
%           gets the undisturbed bathymetry from file file.nc
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
%
if ( nargin ~= 4 | nargout ~= 3 )
  help ksection;
  return
end
%
% Turn off warnings from netCDf:
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
[name,nz]=ncmex('diminq',ncid,'z');
%
% Get 2D or 3D variables:
%
if(itime > 0 & kindex > 0)
 var=ncmex('varget',ncid,vname,[(itime-1) (kindex-1) 0 0],[1 1 -1 -1],1); 
elseif(itime > 0 & kindex < 0)
 var=ncmex('varget',ncid,vname,[(itime-1) 0 0],[1 -1 -1],1);
elseif(itime < 0 & kindex > 0)
 var=ncmex('varget',ncid,vname,[(kindex-1) 0 0],[1 -1 -1],1); 
else
 var=ncmex('varget',ncid,vname,[0 0],[-1 -1],1); 
end
%
if strcmp(vname,'u')
%
% Average velocities to cell centres and extrapolate end value:
%
  var=[ (var(1:nx-1,:)+var(2:nx,:))/2 ; var(nx,:) ];
elseif strcmp(vname,'v')
%
% Average velocities to cell centres and extrapolate end value:
%
  var=[ (var(:,1:ny-1)+var(:,2:ny))/2 var(:,ny) ];
end
%
% Get mask:
%
fsm=ncmex('varget',ncid,'fsm',[0 0],[-1 -1],1);
land=find(fsm==0);
%
var(land)=var(land)*NaN;
%
% Get horizontal coordinates:
%
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
% Close file:
%
ncmex('close',ncid);
