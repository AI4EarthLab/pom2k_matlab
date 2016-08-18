function [art,aru,arv,fsm,dum,dvm]=areas_masks(im,jm,dx,dy,h)

% Set to zero
art=zeros(im,jm);
aru=zeros(im,jm);
arv=zeros(im,jm);
fsm=zeros(im,jm);
dum=zeros(im,jm);
dvm=zeros(im,jm);

% Calculate areas of "t" and "s" cells:
art=dx.*dy;

% Calculate areas of "u" and "v" cells:
% % % for j=2:jm
% % % 	for i=2:im
% % %    		aru(i,j)=.25*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j));
% % %         arv(i,j)=.25*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1));
% % %     end
% % % end
% % % 
% % % for j=1:jm
% % % 	aru(1,j)=aru(2,j);
% % %     arv(1,j)=arv(2,j);
% % % end
% % %       
% % % for i=1:im
% % % 	aru(i,1)=aru(i,2);
% % % 	arv(i,1)=arv(i,2);
% % % end

aru(2:im,2:jm) = .25*(dx(2:im,2:jm)+dx(1:im-1,2:jm)).*(dy(2:im,2:jm)+dy(1:im-1,2:jm));
arv(2:im,2:jm) = .25*(dx(2:im,2:jm)+dx(2:im,1:jm-1)).*(dy(2:im,2:jm)+dy(2:im,1:jm-1));

aru(1,:)=aru(2,:);
aru(:,1)=aru(:,2);

arv(1,:)=arv(2,:);
arv(:,1)=arv(:,2);

%     Initialise and set up free surface mask:
fsm(h>1.0)=1.0;

%     Set up velocity masks:
% % % for j=2:jm
% % % 	for i=2:im
% % %     	dum(i,j)=fsm(i,j)*fsm(i-1,j);
% % %         dvm(i,j)=fsm(i,j)*fsm(i,j-1);
% % %     end
% % % end

dum(2:im,2:jm)=fsm(2:im,2:jm) .* fsm(1:im-1,2:jm);
dvm(2:im,2:jm)=fsm(2:im,2:jm) .* fsm(2:im,1:jm-1);
