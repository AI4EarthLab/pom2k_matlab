function plotall
%
% plotall:         general plotting program, which assembles a set of POMVIZ
%                  functions
%
% Usage: plotall
%
% On harp at the Antarctic CRC of the University of Tasmania, the following
% paths are required:
%
% /opt/matlab/toolbox/matlab/general
% /opt/matlab/toolbox/matlab/ops
% /opt/matlab/toolbox/matlab/lang
% /opt/matlab/toolbox/matlab/elmat
% /opt/matlab/toolbox/matlab/elfun
% /opt/matlab/toolbox/matlab/graph2d
% /opt/matlab/toolbox/matlab/graph3d
% /opt/matlab/toolbox/matlab/specgraph
% /opt/matlab/toolbox/matlab/graphics
% /opt/matlab/toolbox/matlab/uitools
% /opt/matlab/toolbox/matlab/strfun
% /opt/matlab/toolbox/matlab/iofun
% /opt/matlab/toolbox/matlab/datatypes
% /opt/matlab/toolbox/netcdf_matlab5
% /u/johunter/matlab/netcdf
%
% Initial version, JRH 11/12/2001
% Port Latta version, JRH 10/1/2002
% Converted back to general POM version, 21/1/2002
%
% Required plots (0: no plot; 1: plot):
%
profplot = 1;
depthplot=1;
elevplot = 1;
tempplot = 1;
salplot = 1;
%
% Set Run identifier:
%
run='pom2k';
%
% Set file name:
%
file='pom2k.nc';
%
% Set time index:
%
itime=3;
%
% Set i value for section:
%
iindex=32;
%
% Set j value for section:
%
jindex=24;
%
% Set k value for section:
%
kindex=15;
%
% Set height for section:
%
height=-1000;
%
% Set position for velocity key:
%
east_key=5000;
north_key=5000;
%
% Set velocity for key vector:
%
velkey=.2;
%
% Subsample number for vectors:
%
isub=1;
%
% Set scale factor for vectors:
%
scalfac=20000;
%
% Cell indices for profile plots:
%
ij = [ 16   12;
       32   24;
       48   36];
%
%.............................................................................
%
% Profiles:
%
if profplot == 1
%
% Plot series of profiles:
%
  figure;
  getpageprofvel(ij,file,itime);
%
%  print ( '-dpsc2', strcat(run,'_velprof.ps') );
%
% Plot map of cells for profiles:
%
  figure;
  getplot_profloc(ij,file);
%
%  print ( '-dpsc2', strcat(run,'_velprofloc.ps') );
%
end
%
%.............................................................................
%
% Depth:
%
if depthplot==1
%
  figure;
%
  getplot_ksection(file,'h',-1,-1,[0 5000],'Depth (m)')
%
%  print ( '-dpsc2', strcat(run,'_x_y_h.ps') );
%
end
%
%.............................................................................
%
% Elevation:
%
if elevplot==1
%
  figure;
%
  getplot_ksectionvel(file,'elb',itime,-1,itime,kindex,...
    east_key,north_key,velkey,isub,scalfac,[-.5 0.5],'Elevation (m)')
%
%  print ( '-dpsc2', strcat(run,'_x_y_e.ps') );
%
end
%
%.............................................................................
%
% Temperature:
%
if tempplot==1
%
  figure;
%
  getplot_isection(file,'t',itime,iindex,.02,[5 20],'Temperature (deg. C)')
%
%  print ( '-dpsc2', strcat(run,'_y_z_t.ps') );
%
  figure;
%
  getplot_jsection(file,'t',itime,jindex,.02,[5 20],'Temperature (deg. C)')
%
% print ( '-dpsc2', strcat(run,'_x_z_t.ps') );
%
  figure;
%
  getplot_ksectionvel(file,'t',itime,kindex,itime,kindex,...
    east_key,north_key,velkey,isub,scalfac,[5 20],'Temperature (deg. C)')
%
%  print ( '-dpsc2', strcat(run,'_x_y_t_layer.ps') );
%
  figure;
%
  getplot_horiz_sectionvel(file,'t',itime,itime,height,...
    east_key,north_key,velkey,isub,scalfac,[5 20],'Temperature (deg. C)')
%
%  print ( '-dpsc2', strcat(run,'_x_y_t_height.ps') );
%
end
%
%.............................................................................
%
% Salinity:
%
if salplot==1
%
  figure;
%
  getplot_isection(file,'s',itime,iindex,.02,[34 36],'Salinity (PSS)')
%
%  print ( '-dpsc2', strcat(run,'_y_z_s.ps') );
%
  figure;
%
  getplot_jsection(file,'s',itime,jindex,.02,[34 36],'Salinity (PSS)')
%
% print ( '-dpsc2', strcat(run,'_x_z_s.ps') );
%
  figure;
%
  getplot_ksectionvel(file,'s',itime,kindex,itime,kindex,...
    east_key,north_key,velkey,isub,scalfac,[34 36],'Salinity (PSS)')
%
%  print ( '-dpsc2', strcat(run,'_x_y_s_layer.ps') );
%
  figure;
%
  getplot_horiz_sectionvel(file,'s',itime,itime,height,...
    east_key,north_key,velkey,isub,scalfac,[34 36],'Salinity (PSS)')
%
%  print ( '-dpsc2', strcat(run,'_x_y_s_height.ps') );
%
end
%
%.............................................................................
%

