function [var,east_c,north_c]=horiz_section(file,vname,itime,height)
%
% horiz_section:   returns a horizontal section of a 3-D variable at
%                  height, height, from a *.nc file from pom2k
%
% Usage: [var,east_c,north_c]=horiz_section(file,vname,itime,height)
%
% where: var ..... the 2-D section
%        east_c .. horizontal coordinates of corners of model "cell" (metres)
%        north_c . horizontal coordinates of corners of model "cell" (metres)
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%                  (insert negative value if not required)
%        height .. height (up) relative to surface (metres)
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
% "lev=diff(sigma < sigmaheight,1,1)" changed to 
%   "lev=diff(sigma <= sigmaheight,1,1)", JRH 02/01/2002
%
if (nargin ~= 4 |  nargout ~= 3 )
  help horiz_section;
  return
end
%
% Make gross check on height:
%
if height > 10
  help horiz_section;
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
% Get 3-D slab at given time level (if applicable):
if itime < 0
  var3d=ncmex('varget',ncid,vname,[0 0 0],[-1 -1 -1],1);
else
  var3d=ncmex('varget',ncid,vname,[(itime-1) 0 0 0],[1 -1 -1 -1],1);
end
%
% Get mask:
%
fsm=ncmex('varget',ncid,'fsm',[0 0],[-1 -1],1);
land=find(fsm==0);
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
% Get depth and surface elevation:
%
h=ncmex('varget',ncid,'h',[0 0],[-1 -1]);
elb=ncmex('varget',ncid,'elb',[(itime-1) 0 0],[1 -1 -1],1);
%
% Get sigma:
%
if strcmp(vname,'w')
  sigma=ncmex('varget',ncid,'z',0,nz,1);
  sigma=sigma';
  nsigma=nz;
else
  sigma=ncmex('varget',ncid,'zz',0,nz,1);
  sigma=sigma';
  sigma=[ 0 ; sigma(1:nz-1) ; -1 ];
  nsigma=nz+1;
end
%
% Find value of sigma at each point equivalent to depth height:
%
sigmaheight = (height-elb)./(h+elb);
%
% Mask out land:
%
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
lev=diff(sigma <= sigmaheight,1,1);
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
if strcmp(vname,'u')
%
% Average velocities to cell centres and extrapolate end values:
%
  var3d=cat(1,var3d(2,:,:),(var3d(2:nx-1,:,:)+var3d(3:nx,:,:))/2,var3d(nx,:,:));
elseif strcmp(vname,'v')
%
% Average velocities to cell centres and extrapolate end values:
%
  var3d=cat(2,var3d(:,2,:),(var3d(:,2:ny-1,:)+var3d(:,3:ny,:))/2,var3d(:,ny,:));
end
%
% If not w generate top and bottom values:
%
if ~strcmp(vname,'w')
  var3d = cat(3,var3d(:,:,1),var3d(:,:,1:nz-1),var3d(:,:,nz-1));
end
%
% Shift dimension so that leading dimension is vertical:
%
var3d=shiftdim(var3d,2);
%
% Do linear interpolation in sigma:
var=zeros(nx,ny);
%
var(indh)=var3d(ind)+(var3d(ind+1)-var3d(ind)).*(sigmaheight(ind)-sigma(ind))./(sigma(ind+1)-sigma(ind));
%
% Mask:
var(land)=NaN;
var(outside)=NaN;
%
% Close file:
%
ncmex('close',ncid);
