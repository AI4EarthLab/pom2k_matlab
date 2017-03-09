function [u,v,ub,vb]=smoth_update(uf,vf,u,v,ub,vb)
global smoth
ub=u+0.5*smoth*(uf+ub-2.0*u);
u=uf;
vb=v+0.5*smoth*(vf+vb-2.0*v);
v=vf;

% ub=u+0.5*smoth*(ub-2.0*u+uf);
% u=uf;
% vb=v+0.5*smoth*(vb-2.0*v+vf);
% v=vf;
end