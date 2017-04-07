clc;
clear all;
cartesian=0;
problem=input('Input the problem name:','s');
fclim_flag = 0;
grav=9.8060;       % gravity constant (S.I. units)
%Initial velocity:
uvel = 0.2e0;
vvel = 0.e0;
%=============================================================================
[dx,dy,z,zz,cor,rot,h,fsm,dum,dvm,t,s,tclim,sclim,wusurf,wvsurf,wtsurf,  ...
            wssurf,swrad,vflux,e_atmos]=read_basic_variables(problem);
%=============================================================================
im=size(t,1);
jm=size(t,2);
kb=size(zz,2);
imm1=im-1;
jmm1=jm-1;
kbm1=kb-1;
%Initial arrays:
ele=zeros(1,jm)      ;eln=zeros(1,im)      ;els=zeros(1,im)      ;elw=zeros(1,jm)      ;
sbe=zeros(jm,kb)     ;sbn=zeros(im,kb)     ;sbs=zeros(im,kb)     ;sbw=zeros(jm,kb)     ;
tbe=zeros(jm,kb)     ;tbn=zeros(im,kb)     ;tbs=zeros(im,kb)     ;tbw=zeros(jm,kb)     ;
uabe=zeros(1,jm)     ;uabw=zeros(1,jm)     ;ube=zeros(jm,kb)     ;ubw=zeros(jm,kb)     ;
vabn=zeros(1,im)     ;vabs=zeros(1,im)     ;vbn=zeros(im,kb)     ;vbs=zeros(im,kb)     ;
%============================================================================
% Compute east_c and north_c with matrix style. 
% The size of L1 is im*im and the size of R1 is jm*jm.
%  L0=[0 0 0 0 0 0 0]     R0=[0 1 1 1 1]
%     [1 0 0 0 0 0 0]        [0 0 1 1 1] 
%     [1 1 0 0 0 0 0]        [0 0 0 1 1]
%     [1 1 1 0 0 0 0]        [0 0 0 0 1]
%     [1 1 1 1 0 0 0]        [0 0 0 0 0]
%     [1 1 1 1 1 0 0]        
%     [1 1 1 1 1 1 0]        
%tic;
L0=tril(ones(im,im)) - eye(im,im);
R0=triu(ones(jm,jm)) - eye(jm,jm);
east_c=L0*dx;
north_c=dy*R0;

%     The following works properly for any grid:
%
%     Elevation points:

% Define additional matrix
% L =[1  1  0  0  0  0  0  0]    R= [1  0  0  0  0]
%    [0  1  1  0  0  0  0  0]       [1  1  0  0  0]
%    [0  0  1  1  0  0  0  0]       [0  1  1  0  0]
%    [0  0  0  1  1  0  0  0]       [0  0  1  1 -1]
%    [0  0  0  0  1  1  0  0]       [0  0  0  1  3]
%    [0  0  0  0  0  1  1  0]       [0  0  0  0  1]
%    [0  0  0  0  0 -1  3  1]  
%
% T =[0    0    0    0    0   ]
%    [0    0    0    0    0   ]
%    [0    0    0    0    0   ]
%    [0    0    0    0    0   ]
%    [0    0    0    0   -0.5 ]
%    [0    0    0    0    1   ]
%    [0    0   -0.5  1    0   ]
L1=[eye(im) zeros(im,1)];
L2=[zeros(im,1) eye(im)];
L3= zeros(im,im+1); L3(im,im-1)=-1;
L4= zeros(im,im+1); L4(im,im  )= 2;
L=L1+L2+L3+L4;

R1=[eye(jm); zeros(1,jm)];
R2=[zeros(1,jm); eye(jm)];
R3= zeros(jm+1,jm); R3(jm-1,jm)=-1;
R4= zeros(jm+1,jm); R4(jm  ,jm)= 2;
R=R1+R2+R3+R4;

T=zeros(im,jm); T(im,jm-1)=1;  T(im,jm-2)=-0.5;  T(im-1,jm)=1; T(im-2,jm)=-0.5;

% Compute east_e with matrix style
A1=zeros(im+1,jm+1);
A1(1:im,1:jm)=east_c;
B1=0.25*L*A1*R;
B1(im,jm) = sum(sum(T.*B1,1),2);
east_e=B1;

% Compute north_e with matrix style
A2=zeros(im+1,jm+1);
A2(1:im,1:jm)=north_c;
B2=0.25*L*A2*R;
B2(im,jm) = sum(sum(T.*B2,1),2);
north_e=B2;


% Compute east_u with matrix style
A3=zeros(im,jm+1);
A3(1:im,1:jm)=east_c;
B3=A3*R/2.0;
east_u=B3;

% Compute north_u with matrix style
A4=zeros(im,jm+1);
A4(1:im,1:jm)=north_c;
B4=A4*R/2.0;
north_u=B4;


% Compute east_v with matrix style
A5=zeros(im+1,jm);
A5(1:im,1:jm)=east_c;
B5=L*A5/2.0;
east_v=B5;

% Compute north_v with matrix style
A6=zeros(im+1,jm);
A6(1:im,1:jm)=north_c;
B6=L*A6/2.0;
north_v=B6;
%===========================================================================
if cartesian
[z2]=read_cart_z;
[tb,sb]=ztosig(t,s,h,zz,z2);
else
    tb=t;
    sb=s;
end
%===========================================================================
% derived vertical grid variables
dz(1:kb-1)=z(1:kb-1)-z(2:kb);
dzz(1:kb-1)=zz(1:kb-1)-zz(2:kb);
dz(kb)=0.e0;
dzz(kb)=0.e0;
% Compute cell areas 
art=dx.*dy;
aru(2:im,2:jm)=.25e0*(dx(2:im,2:jm)+dx(1:imm1,2:jm)).*(dy(2:im,2:jm)+dy(1:imm1,2:jm));
arv(2:im,2:jm)=.25e0*(dx(2:im,2:jm)+dx(2:im,1:jmm1)).*(dy(2:im,2:jm)+dy(2:im,1:jmm1));
aru(1,:)=aru(2,:);    arv(1,:)=arv(2,:);
aru(:,1)=aru(:,2);    arv(:,1)=arv(:,2);
%===========================================================================
%Initialize velocity in horizontal direction
for k=1:kbm1
    ub(:,:,k)=uvel*dum(:,:);
    vb(:,:,k)=vvel*dvm(:,:);
end
ub(:,:,kb)=0.e0;
vb(:,:,kb)=0.e0;
uab=uvel*dum;
vab=vvel*dvm;
%Initialize surface evolution
elb=-e_atmos;
etb=-e_atmos;
dt=h-e_atmos;
%Initialize boundary conditions
rfe=1.0;  rfw=1.0; rfn=1.0; rfs=1.0;
uabw(2:jmm1)=uab(2,2:jmm1);
uabe(2:jmm1)=uab(imm1,2:jmm1);
vabs(2:imm1)=vab(2:imm1,2);
vabn(2:imm1)=vab(2:imm1,jmm1);
els(2:imm1)=elb(2:imm1,2);
eln(2:imm1)=elb(2:imm1,jmm1);
A7=[zeros(jmm1-1,1) triu(ones(jmm1-1,jmm1-1)) zeros(jmm1-1,1)];
ele= (-cor(imm1,2:jmm1).*uab(imm1,2:jmm1)./grav.*dy(imm1,1:jmm1-1)) * A7;
elw= (-cor(2,2:jmm1).*uab(2,2:jmm1)./grav.*dy(2,1:jmm1-1)) * A7;
ele(2:jmm1)=(ele(2:jmm1)-ele(jmm1/2)).*fsm(im,2:jmm1);
elw(2:jmm1)=(elw(2:jmm1)-elw(jmm1/2)).*fsm(2, 2:jmm1);
ssurf=sb(:,:,1);
tsurf=tb(:,:,1);
if fclim_flag
    tbe(:,1:kbm1) = tclim(im,:,1:kbm1);  tbw(:,1:kbm1) = tclim(1,:,1:kbm1);
    sbe(:,1:kbm1) = sclim(im,:,1:kbm1);  sbw(:,1:kbm1) = sclim(1,:,1:kbm1);
    tbn(:,1:kbm1) = tclim(:,jm,1:kbm1);  tbs(:,1:kbm1) = tclim(:,1,1:kbm1);
    sbn(:,1:kbm1) = sclim(:,jm,1:kbm1);  sbs(:,1:kbm1) = sclim(:,1,1:kbm1);
else
    tbe(:,1:kbm1) = tb(im,:,1:kbm1); tbw(:,1:kbm1) = tb(1,:,1:kbm1);
    sbe(:,1:kbm1) = sb(im,:,1:kbm1); sbw(:,1:kbm1) = sb(1,:,1:kbm1);
    tbn(:,1:kbm1) = tb(:,jm,1:kbm1); tbs(:,1:kbm1) = tb(:,1,1:kbm1);
    sbn(:,1:kbm1) = sb(:,jm,1:kbm1); sbs(:,1:kbm1) = sb(:,1,1:kbm1);
end
%===========================================================================
write_init_variables();