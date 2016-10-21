function getplot_isection(file,vname,itime,iindex,ar,cax,label)
%
% getplot_isection: gets and plots a vertical section of a 3-D variable 
%                  along i=iindex from a *.nc file from pom2k
%
% Usage: getplot_isection(file,vname,itime,iindex,ar,[cax],[label])
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime ... the time index (first output = 1, etc.)
%        iindex .. i-index along which the section is taken
%        ar ...... aspect ratio of plot (< 1 compresses horizontal axis)
%        cax ..... range of data to plot (autoscales if not supplied)
%        label ... legend for colourbar
%
% Initial version, JRH 11/12/2001
% Title added, JRH 02/01/2002
%
if ( nargin < 5 | nargin > 7 )
  help getplot_isection;
  return
end
%
[var,ygrid,zgrid] = isection(file,vname,itime,iindex);
%
tit=sprintf('i = %d, itime = %d',iindex,itime);
%
if nargin == 5
  psection(var,ygrid,zgrid,ar,tit);
elseif nargin == 6
  psection(var,ygrid,zgrid,ar,tit,cax);
else
  psection(var,ygrid,zgrid,ar,tit,cax,label);
end
