function [wt,swr]=read_heat_pnetcdf(iint,iheat)
global problem in_path;
n=sprintf('%04d',fix(iint/iheat));
fnc=[in_path,problem,'.heat.',n,'.nc'];
wt = ncread(fnc,'wtsurf');
swr= ncread(fnc,'swrad');
end