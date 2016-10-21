function getplot_ksection(file,vname,itime,kindex,cax,label)
%
% getplot_ksection: gets and plots a "horizontal" section of a 2-D or 3-D
%                  variable at height, height, from a *.nc file from pom2k
%
% Usage: getplot_ksection(file,vname,itime,kindex,[cax],[label])
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%                  (insert negative value if not required)
%        kindex .. k-index along which the section is taken
%                  (insert negative value if not required)
%        cax ..... range of data to plot (autoscales if not supplied)
%        label ... legend for colourbar
%
% Initial version, JRH 11/12/2001
% Title added, JRH 02/01/2002
%
if ( nargin < 4 | nargin > 6 )
  help getplot_ksection; 
  return
end
%
[var,east_c,north_c]=ksection(file,vname,itime,kindex);
%
tit=sprintf('k = %d, itime = %d',kindex,itime);
%
if nargin == 4
  psection(var,east_c,north_c,1,tit);
elseif nargin == 5
  psection(var,east_c,north_c,1,tit,cax);
else
  psection(var,east_c,north_c,1,tit,cax,label);
end
