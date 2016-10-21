function getplot_ksectionvel(file,vname,itime_1,kindex_1,itime_2,kindex_2,...
  east_key,north_key,velkey,isub,scalfac,cax,label,color)
%
% getplot_ksectionvel: gets and plots a "horizontal" section of a 2-D or
%                  3-D variable, and the horizontal velocity, along k=kindex 
%                  from a *.nc file from pom2k
%
% Usage: getplot_ksectionvel(file,vname,itime_1,kindex_1,itime_2,kindex_2,...
% east_key,north_key,velkey,isub,scalfac,[cax],[label],[color])
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime_1 . the time index (first output = 1, etc.)
%                  (for 2-D or 3-D variable; insert negative value if not required)
%        kindex_1 .. k-index along which the section is taken
%                  (for 2-D or 3-D variable; insert negative value if not required)
%        itime_2 . the time index (first output = 1, etc.) (for velocity)
%        kindex_2 .. k-index along which the section is taken (for velocity)
%        east_key .. east-position of key arrow
%        north_key . north-position of key arrow
%        velkey .. velocity for key arrow (metres/second)
%        isub .... subsample factor
%        scalfac . scale factor for arrows
%        cax ..... range of data to plot (autoscales if not supplied)
%        label ... legend for colourbar
%        color ... colour for arrows
%
% Initial version, JRH 11/12/2001
% Title added, JRH 02/01/2002
%
if ( nargin < 11 | nargin > 14 )
  help getplot_ksection;
  return
end
%
[var,east_c,north_c]=ksection(file,vname,itime_1,kindex_1);
%
tit=sprintf('k = %d, itime = %d',kindex_2,itime_2);
%
if nargin == 11
  psection(var,east_c,north_c,1,tit);
elseif nargin == 12
  psection(var,east_c,north_c,1,tit,cax);
else
  psection(var,east_c,north_c,1,tit,cax,label);
end
%
[vel,east_e,north_e]=ksectionvel(file,itime_2,kindex_2);
%
if nargin < 14
  psectionvel(vel,east_e,north_e,east_key,north_key,velkey,isub,scalfac);
else
  psectionvel(vel,east_e,north_e,east_key,north_key,velkey,isub,scalfac,color);
end
