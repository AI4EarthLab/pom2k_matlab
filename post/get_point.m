function [var]=get_point(file,vname,iindex,jindex,kindex,itime)
%
% get_point:       returns the value of a 2-D or 3-D variable 
%                  at a given point from a *.nc file from pom2k
%
% Usage: [var]=get_point(file,vname,iindex,jindex,kindex,itime)
%
% where: var ..... the required value
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        iindex .. i-index of point
%        jindex .. j-index of point
%        kindex .. k-index of point(increasing downwards from 1 at the surface)
%                  (use negative value if not required)
%        itime ... the time index (first output = 1, etc.)
%                  (use negative value if not required)
%    
% Acknowledgement: this function draws heavily on the "omviz" package
% written by Rich Signell (USGS)
%
% Initial version, JRH 08/01/2002
% Header corrected, JRH 15/01/2002
%
if ( nargin ~= 6 )
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
% Get 2D or 3D variable:
%
if(itime > 0 & kindex > 0)
 var=ncmex('varget',ncid,vname,[(itime-1) (kindex-1) (jindex-1) (iindex-1)],[1 1 1 1],1); 
elseif(itime > 0 & kindex < 0)
 var=ncmex('varget',ncid,vname,[(itime-1) (jindex-1) (iindex-1)],[1 1 1],1);
elseif(itime < 0 & kindex > 0)
 var=ncmex('varget',ncid,vname,[(kindex-1) (jindex-1) (iindex-1)],[1 1 1],1); 
else
 var=ncmex('varget',ncid,vname,[(jindex-1) (iindex-1)],[1 1],1); 
end
%
% Close file:
%
ncmex('close',ncid);
