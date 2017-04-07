function [wu,wv]=read_wind_pnetcdf(iint,iwind)
global problem in_path;
n=sprintf('%04d',fix(iint/iwind));
fnc=[in_path,problem,'.wind.',n,'.nc'];
wu = ncread(fnc,'wusurf')   ;
wv = ncread(fnc,'wvsurf')   ;
end
