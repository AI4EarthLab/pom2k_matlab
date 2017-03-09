function [egb,etb,et,dt,utb,vtb,vfluxb,dt_3d]=internal_update(egf,et,etf,utf,vtf,vfluxf)
global h kb;
        egb=egf;    etb=et;     et=etf;         dt=h+et;
        utb=utf;    vtb=vtf;    vfluxb=vfluxf;  dt_3d=repmat(dt,1,1,kb);
return