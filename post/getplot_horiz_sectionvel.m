function getplot_horiz_sectionvel(file,vname,itime_1,itime_2,height,...
  east_key,north_key,velkey,isub,scalfac,cax,label,color)
%
% getplot_horiz_sectionvel: gets and plots a horizontal section of a 3-D
%                  variable, and the horizontal velocity, at height, height, 
%                  from a *.nc file from pom2k
%
% Usage: getplot_horiz_sectionvel(file,vname,itime_1,itime_2,height,...
%  east_key,north_key,velkey,isub,scalfac,[cax],[label],[color])
%
%        file .... the name of the netCDF file
%        vname ... the name of the netCDF variable
%        itime_1 . the time index (first output = 1, etc.)
%                  (for 3-D variable; insert negative value if not required)
%        itime_2 . the time index (first output = 1, etc.) (for velocity)
%        height .. height (up) relative to surface (metres)
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
if ( nargin < 10 | nargin > 13 )
  help getplot_horiz_sectionvel;
  return
end
%
% Make gross check on height:
%
if height > 10
  help getplot_horiz_sectionvel;
  return
end
%
[var,east_c,north_c]=horiz_section(file,vname,itime_1,height);
%
tit=sprintf('height = %d m, itime = %d',height,itime_2);
%
if nargin == 10
  psection(var,east_c,north_c,1,tit);
elseif nargin == 11
  psection(var,east_c,north_c,1,tit,cax);
else
  psection(var,east_c,north_c,1,tit,cax,label);
end
%
[vel,east_e,north_e]=horiz_sectionvel(file,itime_2,height);
%
if nargin < 13
  psectionvel(vel,east_e,north_e,east_key,north_key,velkey,isub,scalfac);
else
  psectionvel(vel,east_e,north_e,east_key,north_key,velkey,isub,scalfac,color);
end
