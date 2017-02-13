function [art,aru,arv,fsm,dum,dvm]=areas_masks(dx,dy,h)
%load('grid.mat');

global im jm;

% Set to zero
art=zeros(im,jm); aru=zeros(im,jm); arv=zeros(im,jm); 
fsm=zeros(im,jm); dum=zeros(im,jm); dvm=zeros(im,jm);

% Calculate areas of "t" and "s" cells:
art = dx .* dy;
aru = AXB(dx) .* AXB(dy);
arv = AYB(dx) .* AYB(dy);
aru(1,:)=aru(2,:);
aru(:,1)=aru(:,2);
arv(1,:)=arv(2,:);
arv(:,1)=arv(:,2);

% Initialise and set up free surface mask:
fsm(h>1.0)=1.0;

dum(2:im,2:jm)=fsm(2:im,2:jm) .* fsm(1:im-1,2:jm);
dvm(2:im,2:jm)=fsm(2:im,2:jm) .* fsm(2:im,1:jm-1);