function getplot_profvel(file,itime,i,j,iplot,ileg)
%
% getplot_profvel: gets and plots a velocity profile at cell (i,j) from a 
%                  *.nc file from pom2k
%
% Usage: getplot_profvel(file,itime,i,j,iplot,ileg)
%
% where: file .... the name of the netCDF file
%        itime ... the time index (first output = 1, etc.)
%        i ....... i-index of cell
%        j ....... j-index of cell
%        iplot ... subplot number
%        ileg .... 1 for legend, otherwise zero
%
% Initial version, JRH 11/12/2001
% itime added to pprofvel for title, JRH 02/01/2002
%
if ( nargin ~= 6 )
  help getplot_profvel;
  return
end
%
[vel,height,heightlim]=profvel(file,itime,i,j);
%
pprofvel(vel,height,heightlim,i,j,itime,iplot,ileg);
