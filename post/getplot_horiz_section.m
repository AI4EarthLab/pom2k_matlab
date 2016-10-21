function getplot_horiz_section(file,vname,itime,height,cax,label)
%
% getplot_horiz_section: gets and plots a horizontal section of a 3-D
%                  variable at height, height, from a *.nc file from pom2k
%
% Usage: getplot_horiz_section(file,vname,itime,height,[cax],[label])
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%                  (insert negative value if not required)
%        height .. height (up) relative to surface (metres)
%        cax ..... range of data to plot (autoscales if not supplied)
%        label ... legend for colourbar
%
% Initial version, JRH 11/12/2001
% Title added, JRH 02/01/2002
%
if ( nargin < 4 | nargin > 6 )
  help getplot_horiz_section;
  return
end
%
% Make gross check on height:
%
if height > 10
  help getplot_horiz_section;
  return
end
%
[var,east_c,north_c]=horiz_section(file,vname,itime,height);
%
tit=sprintf('height = %d m, itime = %d',height,itime);
%
if nargin == 4
  psection(var,east_c,north_c,1,tit);
elseif nargin == 5
  psection(var,east_c,north_c,1,tit,cax);
else
  psection(var,east_c,north_c,1,tit,cax,label);
end
