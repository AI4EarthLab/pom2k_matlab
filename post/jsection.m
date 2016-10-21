function [var,xgrid,zgrid] = jsection(file,vname,itime,jindex)
%
% jsection:        returns a vertical section of a 3-D variable along
%                  j=jindex from a *.nc file from pom2k
%
% Usage: [var,xgrid,zgrid] = jsection(file,vname,itime,jindex)
%
% where: var ..... the 2-D section
%        xgrid ... horizontal coordinates of corners of model "cell", measured
%                  along a line of constant j (metres)
%        zgrid ... vertical coordinates of corners of model "cell" (metres)
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%        jindex .. j-index along which the section is taken
%
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 11/12/2001
%
if (nargin ~= 4 | nargout ~= 3)
  help jsection;
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
  if jindex == ny
    var=ncmex('varget',ncid,vname,[(itime-1) 0 (jindex-1) 0],...
             [1 nz-1 1 nx],1);
  else
%
% Average to cell centers:
%
    var(:,:,1)=ncmex('varget',ncid,vname,[(itime-1) 0 (jindex-1) 0],...
                    [1 nz-1 1 nx],1);
    var(:,:,2)=ncmex('varget',ncid,vname,[(itime-1) 0 jindex 0],...
                    [1 nz-1 1 nx],1);
    var=sum(var,3)/2;
  end
else
  var=ncmex('varget',ncid,vname,[(itime-1) 0 (jindex-1) 0],...
           [1 nz-1 1 nx],1);
end
var=reshape(var,nx,nz-1);
%
h=ncmex('varget',ncid,'h',[(jindex-1) 0],[1 nx],1);
%
elb=ncmex('varget',ncid,'elb',[(itime-1) (jindex-1) 0],[1 1 nx],1);
%
fsm=ncmex('varget',ncid,'fsm',[(jindex-1) 0],[1 nx],1);
%
ind=find(fsm==0);
h(ind)=h(ind)*0;
%
if strcmp(vname,'v')
%
% Average velocities to cell centres and extrapolate end value:
%
  var=[ (var(1:nx-1,:)+var(2:nx,:))/2 ; var(nx,:) ];
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
dx=ncmex('varget',ncid,'dx',[(jindex-1) 0],[1 nx],1);
%
% Integrate dx:
%
xgrid=[0 ; cumsum(dx,1)];
%
if h(1) == 0
  hstart=0;
elseif h(2) == 0
  hstart=h(1);
else
  hstart=1.5*h(1)-0.5*h(2);
end
%
if h(nx) == 0
  hstop=0;
elseif h(nx-1) == 0
  hstop=h(nx);
else
  hstop=1.5*h(nx)-0.5*h(nx-1);
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
if elb(nx) == 0
  elbstop=0;
elseif elb(nx-1) == 0
  elbstop=elb(nx);
else
  elbstop=1.5*elb(nx)-0.5*elb(nx-1);
end
%
% Don't do average when one depth is zero:
%
for i = 1:nx-1
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
h=h(1:nx-1);
elb=elb(1:nx-1);
h=[ hstart ; h ; hstop ];
elb=[ elbstart ; elb ; elbstop ];
%
[m n]=size(sigma);
xgrid=xgrid*ones(1,m);
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
