function [vel,east_e,north_e]=horiz_sectionvel(file,itime,height)
%
% horiz_sectionvel:returns a horizontal section of horizontal velocity
%                  vectors at height, height, from a *.nc file from pom2k
%
% Usage: [vel,east_e,north_e]=horiz_sectionvel(file,itime,height)
%
% where: vel ..... complex 2-D section of horizontal velocity vectors, with 
%                  components in east and north
%        east_e .. horizontal coordinates of centres of model "cell" (metres)
%        north_e . horizontal coordinates of centres of model "cell" (metres)
%
%        file .... the name of the netCDF file
%        itime ... the time index (first output = 1, etc.)
%        height .. height (up) relative to surface (metres)
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
% Sign of rot corrected, JRH 19/12/2001
%
if (nargin ~= 3 | nargout ~= 3),
  help horiz_sectionvel;
  return
end
%
% Make gross check on height:
%
if height > 10
  help horiz_sectionvel;
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
% Get 3-D slab of velocity at given time level:
%
u=ncmex('varget',ncid,'u',[(itime-1) 0 0 0],[1 -1 -1 -1],1);
v=ncmex('varget',ncid,'v',[(itime-1) 0 0 0],[1 -1 -1 -1],1);
%
% Get mask:
%
fsm=ncmex('varget',ncid,'fsm',[0 0],[-1 -1],1);
land=find(fsm==0);
%
% Get horizontal coordinates and rotation angle:
%
east_e=ncmex('varget',ncid,'east_e',[0 0],[-1 -1],1);
north_e=ncmex('varget',ncid,'north_e',[0 0],[-1 -1],1);
rot=ncmex('varget',ncid,'rot',[0 0],[-1 -1],1);
%
% Get depth and surface elevation:
%
h=ncmex('varget',ncid,'h',[0 0],[-1 -1]);
elb=ncmex('varget',ncid,'elb',[(itime-1) 0 0],[1 -1 -1],1);
%
% Get sigma:
%
sigma=ncmex('varget',ncid,'zz',0,nz,1);
sigma=sigma';
sigma=[ 0 ; sigma(1:nz-1) ; -1 ];
nsigma=nz+1;
%
% Find value of sigma at each point equivalent to depth height:
%
sigmaheight = (height-elb)./(h+elb);
%
% Mask out land:
sigmaheight(land)=NaN;
%
% Make mask for values outside water column:
%
outside=find( sigmaheight < -1 | sigmaheight > 0 );
sigmaheight(outside)=NaN;
%
% Generate 3-D "sigma" arrays, with leading dimension vertical:
%
sigma=reshape(sigma*ones(1,nx*ny),nsigma,nx,ny);
%
sigmaheight=reshape(ones(nsigma,1)*sigmaheight(:).',nsigma,nx,ny);
%
lev=diff(sigma < sigmaheight,1,1);
lev(nsigma,:,:)=reshape(zeros(1,nx*ny),1,nx,ny);
%
% Indices of sigma levels immediately above sigmaheight:
%
ind=find(lev == 1);
%
% Horizontal indices where interpolation required:
%
indh=find(sum(lev,1));
%
% Average velocities to cell centres and extrapolate end values:
%
u=cat(1,u(2,:,:),(u(2:nx-1,:,:)+u(3:nx,:,:))/2,u(nx,:,:));
v=cat(2,v(:,2,:),(v(:,2:ny-1,:)+v(:,3:ny,:))/2,v(:,ny,:));
%
% Generate complex horizontal velocity in (x,y,z):
%
uv3d=u+sqrt(-1)*v;
%
% Generate top and bottom values:
uv3d = cat(3,uv3d(:,:,1),uv3d(:,:,1:nz-1),uv3d(:,:,nz-1));
%
% Shift dimension so that leading dimension is vertical:
%
uv3d=shiftdim(uv3d,2);
%
% Do linear interpolation in sigma:
%
uv=zeros(nx,ny);
%
uv(indh)=uv3d(ind)+(uv3d(ind+1)-uv3d(ind)).*(sigmaheight(ind)-sigma(ind))./(sigma(ind+1)-sigma(ind));
%
% Mask:
%
uv(land)=NaN;
uv(outside)=NaN;
%
% Rotate into east/north components:
%
vel=uv.*exp(sqrt(-1)*rot);
%
% Close file:
%
ncmex('close',ncid);
