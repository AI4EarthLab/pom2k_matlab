function [art,aru,arv,fsm,dum,dvm]=areas_masks(art,aru,arv,fsm,dum,dvm,...
        im,jm,dx,dy,h)


%     Calculate areas of "t" and "s" cells:
art=dx.*dy;
%     Calculate areas of "u" and "v" cells:


for j=2:jm
	for i=2:im
   		aru(i,j)=.25e0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j));
        arv(i,j)=.25e0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1));
    end
end

for j=1:jm
	aru(1,j)=aru(2,j);
    arv(1,j)=arv(2,j);
end
      
for i=1:im
	aru(i,1)=aru(i,2);
	arv(i,1)=arv(i,2);
end


%     Initialise and set up free surface mask:
fsm=zeros(im,jm);
dum=zeros(im,jm);
dvm=zeros(im,jm);
fsm(h>1.e0)=1.e0;

%     Set up velocity masks:
for j=2:jm
	for i=2:im
    	dum(i,j)=fsm(i,j)*fsm(i-1,j);
        dvm(i,j)=fsm(i,j)*fsm(i,j-1);
    end
end

return 
end

