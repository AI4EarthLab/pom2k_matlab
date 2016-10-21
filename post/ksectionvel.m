function [vel,east_e,north_e]=ksectionvel(file,itime,kindex)
%
% ksectionvel:     returns a "horizontal" section of horizontal velocity
%                  vectors along k=kindex from a *.nc file from pom2k
%
% Usage: [vel,east_e,north_e]=ksectionvel(file,itime,kindex)
%
% where: vel ..... complex 2-D section of horizontal velocity vectors, with 
%                  components in east and north
%        east_e .. horizontal coordinates of centres of model "cell" (metres)
%        north_e . horizontal coordinates of centres of model "cell" (metres)
%
%        file .... the name of the netCDF file
%        itime ... the time index (first output = 1, etc.)
%        kindex .. k-index across which the section is taken (increasing
%                  downwards from 1 at the surface)
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
% Sign of rot corrected, JRH 19/12/2001
%
if (nargin ~= 3 | nargout ~= 3),
  help ksectionvel;
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
if(ncid==-1),
  disp(['File ' file ' not found'])
  return
end
%
[name,nx]=ncmex('diminq',ncid,'x');
[name,ny]=ncmex('diminq',ncid,'y');
[name,nz]=ncmex('diminq',ncid,'z');
%
% Get velocity:
%
u=ncmex('varget',ncid,'u',[(itime-1) (kindex-1) 0 0],[1 1 -1 -1],1);
v=ncmex('varget',ncid,'v',[(itime-1) (kindex-1) 0 0],[1 1 -1 -1],1);
%
% Average velocity to cell centres and extrapolate end values:
%
u(2:nx-1,:)=(u(2:nx-1,:)+u(3:nx,:))/2;
u(1,:)=u(2,:);
u(nx,:)=u(nx-1,:);
%
v(:,2:ny-1)=(v(:,2:ny-1)+v(:,3:ny))/2;
v(:,1)=v(:,2);
v(:,ny)=v(:,ny-1);
%
% Generate complex velocity:
%
uv=u+sqrt(-1)*v;
%
% Get mask:
%
fsm=ncmex('varget',ncid,'fsm',[0 0],[-1 -1],1);
land=find(fsm==0);
uv(land)=NaN;
%
% Rotate into east/north components:
%
rot=ncmex('varget',ncid,'rot',[0 0],[-1 -1],1);
vel=uv.*exp(sqrt(-1)*rot);
%
% Get horizontal coordinates:
%
east_e=ncmex('varget',ncid,'east_e',[0 0],[-1 -1],1);
north_e=ncmex('varget',ncid,'north_e',[0 0],[-1 -1],1);
%
% Close file:
%
ncmex('close',ncid);
