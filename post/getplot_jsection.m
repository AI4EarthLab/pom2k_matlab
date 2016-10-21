function getplot_jsection(file,vname,itime,jindex,ar,cax,label)
%
% getplot_jsection: gets and plots a vertical section of a 3-D variable 
%                  along j=jindex from a *.nc file from pom2k
%
% Usage: getplot_jsection(file,vname,itime,jindex,ar,[cax],[label])
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%        jindex .. j-index along which the section is taken
%        ar ...... aspect ratio of plot (< 1 compresses horizontal axis)
%        cax ..... range of data to plot (autoscales if not supplied)
%        label ... legend for colourbar
%
% Initial version, JRH 11/12/2001
% Title added, JRH 02/01/2002
%
if ( nargin < 5 | nargin > 7 )
  help getplot_jsection;
  return
end
%
[var,xgrid,zgrid] = jsection(file,vname,itime,jindex);
%
tit=sprintf('j = %d, itime = %d',jindex,itime);
%
if nargin == 5
  psection(var,xgrid,zgrid,ar,tit);
elseif nargin == 6
  psection(var,xgrid,zgrid,ar,tit,cax);
else
  psection(var,xgrid,zgrid,ar,tit,cax,label);
end
