function[tb,sb,tclim,sclim,ub,uab,elb,etb,dt,aam2d,rho,rmean,wusurf,wvsurf,dt_3d] = seamount(e_atmos,aam)
    
global im jm kb pi slmax kbm1 zz tbias sbias jmm1 imm1 grav dx dy cor h uabe uabw      ...
       rfe rfw ele elw rfn rfs tbe tbw tbn tbs sbe sbw sbn sbs east_c north_c east_e   ...
       north_e east_u north_u east_v north_v rot art aru arv fsm dum dvm dx_3d dy_3d   ...
       cor_3d h_3d art_3d aru_3d arv_3d fsm_3d dum_3d dvm_3d;
       
tb=zeros(im,jm,kb)   ;sb=zeros(im,jm,kb)   ;tclim=zeros(im,jm,kb);sclim=zeros(im,jm,kb);
ub=zeros(im,jm,kb)   ;rmean=zeros(im,jm,kb);wusurf=zeros(im,jm)  ;wvsurf=zeros(im,jm)  ;

% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Sets up for seamount problem.                       *
% *                                                                    *
% **********************************************************************
%
%     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
delh=0.9e0;
%     Grid size:
delx=8000.e0;
%     Radius island or seamount:
ra=25000.e0;
%     Current velocity:
vel=0.2e0;

dx=repmat( (delx-delx*sin(pi*[1:im]'/im)/2.0), 1, jm );
dy=repmat( (delx-delx*sin(pi*[1:jm]/jm)/2.0), im, 1  );
cor(:,:)=1.e-4;

%     Calculate horizontal coordinates of grid points and rotation
%     angle.
%
%     NOTE that this is introduced solely for the benefit of any post-
%     processing software, and in order to conform with the requirements
%     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
%
%     There are four horizontal coordinate systems, denoted by the
%     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
%     "e" is an elevation point and "c" is a cell corner), as shown
%     below. In addition, "east_*" is an easting and "north_*" is a
%     northing. Hence the coordinates of the "u" points are given by
%     (east_u,north_u).
%
%     Also, if the centre point of the cell shown below is at
%     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
%     the coordinates of the western of the two "u" points,
%     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
%     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
%     coordinates of the southwestern corner point of the cell. The
%     southwest corner of the entire grid is at
%     (east_c(1,1),north_c(1,1)).
%
%                      |              |
%                    --c------v-------c--
%                      |              |
%                      |              |
%                      |              |
%                      |              |
%                      u      e       u
%                      |              |
%                      |              |
%                      |              |
%                      |              |
%                    --c------v-------c--
%                      |              |
%
%
%     NOTE that the following calculation of east_c and north_c only
%     works properly for a rectangular grid with east and north aligned
%     with i and j, respectively:
%

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

%     rot is the angle (radians, anticlockwise) of the i-axis relative
%     to east, averaged to a cell centre:
%
%     (NOTE that the following calculation of rot only works properly
%     for this particular rectangular grid)
%
rot=zeros(im,jm);

%     Define depth:

W1= east_c  - repmat( east_c ((im+1)/2,:), im, 1 );
W2= north_c - repmat( north_c(:,(jm+1)/2), 1 , jm);
h=4500.0*(1.e0-delh*exp((-W1.^2 - W2.^2)/ra^2));
h(h<1.0) = 1.0;

%     Close the north and south boundaries to form a channel:
h(:,1) =1.0;
h(:,jm)=1.0;

%     Calculate areas and masks:
%[art,aru,arv,fsm,dum,dvm]=areas_masks(im,jm,dx,dy,h);
areas_masks();
    for j=2:jmm1
        for i=2:imm1
                wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))  ...
                    *0.25*(dvm(i,j+1)+dvm(i-1,j+1)     ...
                    +dvm(i-1,j)+dvm(i,j));
                wvsurf(i,j)=0.e0;
        end
     end
dx_3d=repmat(dx,1,1,kb);
dy_3d=repmat(dy,1,1,kb);
cor_3d=repmat(cor,1,1,kb);
h_3d=repmat(h,1,1,kb);
art_3d=repmat(art,1,1,kb);
aru_3d=repmat(aru,1,1,kb);
arv_3d=repmat(arv,1,1,kb);
fsm_3d=repmat(fsm,1,1,kb);
dum_3d=repmat(dum,1,1,kb);
dvm_3d=repmat(dvm,1,1,kb);


%     Adjust bottom topography so that cell to cell variations
%     in h for not exceed parameter slmax:
if(slmax < 1.e0) 
	[h] = slpmax(h,im,jm,fsm,slmax);
end

% Set initial conditions:
for k=1:kbm1
    tb(:,:,k)=5.0+15.0*exp(zz(k)*h(:,:)/1000.0)-tbias;
    sb(:,:,k)=35.0-sbias;
    tclim(:,:,k)=tb(:,:,k);
    sclim(:,:,k)=sb(:,:,k);
    ub(:,:,k)=vel*dum(:,:);
end 

%     Initialise uab and vab as necessary
%     (NOTE that these have already been initialised to zero in the
%     main program):
uab=vel*dum;

%     Set surface boundary conditions, e_atmos, vflux, wusurf,
%     wvsurf, wtsurf, wssurf and swrad, as necessary
%     (NOTE:
%      1. These have all been initialised to zero in the main program.
%      2. The temperature and salinity of inflowing water must be
%         defined relative to tbias and sbias.):

%     Initialise elb, etb, dt and aam2d:
elb=-e_atmos;
etb=-e_atmos;
dt=h-e_atmos;
aam2d=aam(:,:,1);

dt_3d=repmat(dt,1,1,kb);

 [rho]=dens(sb,tb,h_3d,fsm_3d);

%     Generated horizontally averaged density field (in this
%     application, the initial condition for density is a function
%     of z (the vertical cartesian coordinate) -- when this is not
%     so, make sure that rmean has been area averaged BEFORE transfer
%     to sigma coordinates):

for k=1:kbm1
    rmean(:,:,k) = rho(:,:,k);
end

%     Set lateral boundary conditions, for use in subroutine bcond
%     (in the seamount problem, the east and west boundaries are open,
%     while the south and north boundaries are closed through the
%     specification of the masks fsm, dum and dvm):
rfe=1.e0;
rfw=1.e0;
rfn=1.e0;
rfs=1.e0;

uabw(2:jmm1)=uab(2,2:jmm1);
uabe(2:jmm1)=uab(imm1,2:jmm1);
% Compute ele and elw with matrix style. 
% The size of A7 is (jmm1-1)*jm.
%  A7=[0 1 1 1 0]
%     [0 0 1 1 0] 
%     [0 0 0 1 0]
A7=[zeros(jmm1-1,1) triu(ones(jmm1-1,jmm1-1)) zeros(jmm1-1,1)];
ele= (-cor(imm1,2:jmm1).*uab(imm1,2:jmm1)./grav.*dy(imm1,1:jmm1-1)) * A7;
elw= (-cor(2,2:jmm1).*uab(2,2:jmm1)./grav.*dy(2,1:jmm1-1)) * A7;

%     Adjust boundary elevations so that they are zero in the middle
%     of the channel:
ele(2:jmm1)=(ele(2:jmm1)-ele(jmm1/2)).*fsm(im,2:jmm1);
elw(2:jmm1)=(elw(2:jmm1)-elw(jmm1/2)).*fsm(2, 2:jmm1);

tbe(:,1:kbm1) = tb(im,:,1:kbm1);
tbw(:,1:kbm1) = tb(1,:,1:kbm1);
sbe(:,1:kbm1) = sb(im,:,1:kbm1);
sbw(:,1:kbm1) = sb(1,:,1:kbm1);
tbn(:,1:kbm1) = tb(:,jm,1:kbm1);
tbs(:,1:kbm1) = tb(:,1,1:kbm1);
sbn(:,1:kbm1) = sb(:,jm,1:kbm1);
sbs(:,1:kbm1) = sb(:,1,1:kbm1);
