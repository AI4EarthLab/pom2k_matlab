function [var,ygrid,zgrid] = isection(file,vname,itime,iindex)
%
% isection:        returns a vertical section of a 3-D variable along
%                  i=iindex from a *.nc file from pom2k
%
% Usage: [var,ygrid,zgrid] = isection(file,vname,itime,iindex)
%
% where: var ..... the 2-D section
%        ygrid ... horizontal coordinates of corners of model "cell", measured
%                  along a line of constant i (metres)
%        zgrid ... vertical coordinates of corners of model "cell" (metres)
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%        iindex .. i-index along which the section is taken
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
%
if (nargin ~= 4 | nargout ~= 3)
  help isection;
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
[name,nz]=ncmex('diminq',ncid,'z');
%
if strcmp(vname,'u')
  if iindex == nx
    var=ncmex('varget',ncid,vname,[(itime-1) 0 0 (iindex-1)],...
             [1 nz-1 ny 1],1);
  else
%
% Average to cell centers:
%
    var(:,:,1)=ncmex('varget',ncid,vname,[(itime-1) 0 0 (iindex-1)],...
                    [1 nz-1 ny 1],1);
    var(:,:,2)=ncmex('varget',ncid,vname,[(itime-1) 0 0 iindex],...
                    [1 nz-1 ny 1],1);
    var=sum(var,3)/2;
  end
else
  var=ncmex('varget',ncid,vname,[(itime-1) 0 0 (iindex-1)],...
           [1 nz-1 ny 1],1);
end
var=reshape(var,ny,nz-1);
%
h=ncmex('varget',ncid,'h',[0 (iindex-1)],[ny 1],1);
h=h';
%
elb=ncmex('varget',ncid,'elb',[(itime-1) 0 (iindex-1)],[1 ny 1],1);
elb=elb';
%
fsm=ncmex('varget',ncid,'fsm',[0 (iindex-1)],[ny 1],1);
fsm=fsm';
%
ind=find(fsm==0);
h(ind)=h(ind)*0;
%
if strcmp(vname,'v')
%
% Average velocities to cell centres and extrapolate end value:
%
  var=[ (var(1:ny-1,:)+var(2:ny,:))/2 ; var(ny,:) ];
end
var(ind,:)=var(ind,:)*NaN;
%
% Get sigma:
%
if strcmp(vname,'w')
  sigma=ncmex('varget',ncid,'zz',0,nz,1);
  sigma=sigma';
%
% Add extra line of zeros at bottom 
% (NOTE use of u at surface so as to include NaN where appropriate):
%
  var=[var var(:,1)*0.];
else
  sigma=ncmex('varget',ncid,'z',0,nz,1);
  sigma=sigma';
end
%
% Get grid interval:
%
dy=ncmex('varget',ncid,'dy',[0 (iindex-1)],[ny 1],1);
dy=dy';
%
% Integrate dy:
%
ygrid=[0 ; cumsum(dy,1)];
%
if h(1) == 0
  hstart=0;
elseif h(2) == 0
  hstart=h(1);
else
  hstart=1.5*h(1)-0.5*h(2);
end
%
if h(ny) == 0
  hstop=0;
elseif h(ny-1) == 0
  hstop=h(ny);
else
  hstop=1.5*h(ny)-0.5*h(ny-1);
end
%
if elb(1) == 0
  elbstart=0;
elseif elb(2) == 0
  elbstart=elb(1);
else
  elbstart=1.5*elb(1)-0.5*elb(2);
end
%
if elb(ny) == 0
  elbstop=0;
elseif elb(ny-1) == 0
  elbstop=elb(ny);
else
  elbstop=1.5*elb(ny)-0.5*elb(ny-1);
end
%
% Don't do average when one depth is zero:
%
for i = 1:ny-1
  if h(i) == 0 & h(i+1) == 0
    h(i) = 0;
    elb(i) = 0;
  elseif h(i) == 0
    h(i) = h(i+1);
    elb(i) = elb(i+1);
  elseif h(i+1) == 0
    h(i) = h(i);
    elb(i) = elb(i);
  else
    h(i) = (h(i)+h(i+1))/2;
    elb(i) = (elb(i)+elb(i+1))/2;
  end
end
%
h=h(1:ny-1);
elb=elb(1:ny-1);
h=[ hstart ; h ; hstop ];
elb=[ elbstart ; elb ; elbstop ];
%
[m n]=size(sigma);
ygrid=ygrid*ones(1,m);
%
% Calculate depth:
%
zgrid=(h+elb)*(sigma')+elb*ones(size(sigma'));
%
% Zero any NaN:
%
ind=find(isnan(zgrid));
zgrid(ind)=zeros(size(ind));
%
% Close file:
%
ncmex('close',ncid);
