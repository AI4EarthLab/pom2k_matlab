function [ws]=read_water_pnetcdf(iint,iwater)
global problem in_path;
n=sprintf('%04d',fix(iint/iwater));
fnc=[in_path,problem,'.water.',n,'.nc'];
ws = ncread(fnc,'wssurf');
end