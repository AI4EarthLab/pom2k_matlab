function [art,aru,arv,fsm,dum,dvm]=areas_masks(art,aru,arv,fsm,dum,dvm,im,jm,dx,dy,h)
%     Calculate areas of "t" and "s" cells:
art=dx.*dy;

%     Calculate areas of "u" and "v" cells:
for j=2:jm
	for i=2:im
   		aru(i,j)=.25*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j));
        arv(i,j)=.25*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1));
    end
end

W1=zeros(im,jm);
W2=zeros(im,jm);
W3=zeros(im,jm);
W4=zeros(im,jm);
for j=2:jm
	for i=2:im
        W1(i,j)=(dx(i,j)+dx(i-1,j));
        W2(i,j)=(dy(i,j)+dy(i-1,j));
        W3(i,j)=(dx(i,j)+dx(i,j-1));
        W4(i,j)=(dy(i,j)+dy(i,j-1));  		
    end
end

aru1=0.25*W1.*W2;
arv1=0.25*W3.*W4;

L1=[eye(im-1) zeros(im-1,1)];
L2=[zeros(im-1,1) eye(im-1)];
L=L1+L2;

R1=[eye(jm-1); zeros(1,jm-1)];
R2=[zeros(1,jm-1); eye(jm-1)];
R=R1+R2;

A=0.25*(L*dx).*(L*dy);
B=0.25*(dx*R).*(dy*R);


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
fsm(h>1.0)=1.0;

%     Set up velocity masks:
for j=2:jm
	for i=2:im
    	dum(i,j)=fsm(i,j)*fsm(i-1,j);
        dvm(i,j)=fsm(i,j)*fsm(i,j-1);
    end
end

return 
end

