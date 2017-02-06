function [dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,rot,tatm,satm,vfluxf,tbias,sbias] = box(dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,...
        e_atmos,aam,im,jm,kb,imm1,jmm1,kbm1,slmax,zz,tbias,sbias,grav,rhoref,rot,tatm,satm,vfluxf)
% **********************************************************************
% *                                                                    *
% * FUN%TION    :  Sets up conservation box problem.                   *
% *                                                                    *
% *                This basin uses the same grid as the seamount       *
% *                problem, but it has a flat bottom, is surrounded by *
% *                walls and is initialised with uniform salinity and  *
% *                temperature. It is forced by a surface input of     *
% *                water of the same temperature and salinity as the   *
% *                water in the basin. Therefore, the temperature and  *
% *                salinity in the basin should not change, and the    *
% *                free surface should fall at a rate vflux. It is also*
% *                forced by a steady atmospheric pressure field which *
% *                depresses the southwestern half of the model by 1 m *
% *                and elevates the northeastern half of the model by  *
% *                1 m.                                                *
% *                                                                    *
% *                Since this problem defines its own fixed e_atmos,   *
% *                tatm, satm and e_atmos, comment out corresponding   *
% *                declarations after the for 9000 statement in main    *
% *                program.                                            *
% **********************************************************************


%     Water depth:
      depth=4500.e0;
%     Grid size:
      delx=8000.e0;
%
%     Set up grid dimensions, areas of free surface cells, and
%     %oriolis parameter:
%



for j=1:jm
    for i=1:im
%     For constant grid size:
%         dx(i,j)=delx
%         dy(i,j)=delx
%     Set up grid dimensions, areas of free surface cells, and
%     Coriolis parameter:
%
%     For variable grid size:
          dx(i,j)=delx-delx*sin(pi*i/im)/2.e0;
          dy(i,j)=delx-delx*sin(pi*j/jm)/2.e0;
          cor(i,j)=1.e-4;
	end
end

%     Calculate horizontal coordinates of grid points and rotation
%     angle.
%
%     NOTE that this is introduced solely for the benefit of any post-
%     processing software, and in order to conform with the requirements
%     of the Net%DF %limate and Forecast (%F) Metadata %onventions.
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
for j=1:jm
	east_c(1,j)=0.e0;
	for i=2:im
		east_c(i,j)=east_c(i-1,j)+dx(i-1,j);
		end
	end

	for i=1:im
		north_c(i,1)=0.e0;
		for j=2:jm
			north_c(i,j)=north_c(i,j-1)+dy(i,j-1);
		end
	end
%
%     The following works properly for any grid:
%
%     Elevation points:
%
for j=1:jm-1
	for i=1:im-1
		east_e(i,j)=(east_c(i,j)+east_c(i+1,j)+east_c(i,j+1)+east_c(i+1,j+1))/4.e0;
		north_e(i,j)=(north_c(i,j)+north_c(i+1,j)+north_c(i,j+1)+north_c(i+1,j+1))/4.e0;
	end
end
%
%     Extrapolate ends:
%
for i=1:im-1
	east_e(i,jm)=((east_c(i,jm)+east_c(i+1,jm))*3.e0-east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0;
	north_e(i,jm)=((north_c(i,jm)+north_c(i+1,jm))*3.e0-north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0;
end
%
for j=1:jm-1
	east_e(im,j)=((east_c(im,j)+east_c(im,j+1))*3.e0-east_c(im-1,j)-east_c(im-1,j+1))/4.e0;
	north_e(im,j)=((north_c(im,j)+north_c(im,j+1))*3.e0-north_c(im-1,j)-north_c(im-1,j+1))/4.e0;
end
%
east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)-(east_e(im-2,jm)+east_e(im,jm-2))/2.e0;
north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)-(north_e(im-2,jm)+north_e(im,jm-2))/2.e0;

%     u-points:
for j=1:jm-1
	for i=1:im
		east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0;
		north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0;
	end
end
%     Extrapolate ends:
for i=1:im
	east_u(i,jm)=(east_c(i,jm)*3.0-east_c(i,jm-1))/2.e0;
	north_u(i,jm)=(north_c(i,jm)*3.0-north_c(i,jm-1))/2.e0;
end
 
%     v-points:
for j=1:jm
	for i=1:im-1
		east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0;
		north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0;
	end
end
%     Extrapolate ends:
for j=1:jm
	east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0;
	north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0;
end
 
%     rot is the angle (radians, anticlockwise) of the i-axis relative
%     to east, averaged to a cell centre:
%
%     (NOTE that the following calculation of rot only works properly
%     for this particular rectangular grid)
%
rot=zeros(im,jm);
 
%     Define depth:
for i=1:im
	for j=1:jm
    	h(i,j)=depth;
    end
end
%
%     Close the north and south boundaries:
%
for i=1:im
	h(i,1)=1.e0;
    h(i,jm)=1.e0;
end
 
%     Close the east and west boundaries:
for j=1:jm
	h(1,j)=1.e0;
	h(im,j)=1.e0;
end
 
%     Calculate areas and masks:
[art,aru,arv,fsm,dum,dvm]=areas_masks(art,aru,arv,fsm,dum,dvm,...
        im,jm,dx,dy,h);
 
%     Adjust bottom topography so that cell to cell variations
%     in h for not exceed parameter slmax:
if(slmax < 1.0) 
	[h]=slpmax(h,im,jm,fsm,slmax);
end
 
%     Set tbias and sbias here for test (tbias and sbias would
%     normally only be set in the main program):
tbias=10.e0;
sbias=20.e0;

 
%     Set initial conditions:


for k=1:kbm1
	for j=1:jm
		for i=1:im
			tb(i,j,k)=20.e0-tbias;
			sb(i,j,k)=35.e0-sbias;
			tclim(i,j,k)=tb(i,j,k);
			sclim(i,j,k)=sb(i,j,k);
		end
	end
end
 
%     Initialise uab and vab as necessary
%     (NOTE that these have already been initialised to zero in the
%     main program):


%     Set surface boundary conditions, e_atmos, vflux, wusurf,
%     wvsurf, wtsurf, wssurf and swrad, as necessary
%     (NOTE:
%      1. These have all been initialised to zero in the main program.
%      2. The temperature and salinity of inflowing water must be
%         defined relative to tbias and sbias.):
%
for j=1:jm
	for i=1:im
		if(i+j-57 <= 0)
			e_atmos(i,j)=1.e0;
		else
			e_atmos(i,j)=-1.e0;
		end
%
%     Ensure atmospheric pressure cannot make water depth go negative:
%
		e_atmos(i,j)=min(e_atmos(i,j),h(i,j));
		vfluxf(i,j)=-0.0001e0;
%     See main program, just after "Begin numerical integration", for
%     an explanation of these terms:
		tatm=20.e0;
		satm=35.e0;
	end
end

%     Initialise elb, etb, dt and aam2d:
%
elb=-e_atmos;
etb=-e_atmos;
dt=h-e_atmos;
aam2d=aam(:,:,1);

[sb,tb,rho]=dens(sb,tb,rho,...
    im,jm,kbm1,tbias,sbias,grav,rhoref,zz,h,fsm);
%
%     Generated horizontally averaged density field (in this
%     application, the initial condition for density is a function
%     of z (the vertical cartesian coordinate) -- when this is not
%     so, make sure that rmean has been area averaged BEFORE transfer
%     to sigma coordinates):

for k=1:kbm1
   for j=1:jm
       for i=1:im
          rmean(i,j,k)=rho(i,j,k); 
       end
   end
end
%     Set lateral boundary conditions, for use in subroutine bcond
%     (in this problem, all lateral boundaries are closed through
%     the specification of the masks fsm, dum and dvm):
rfe=1.e0;
rfw=1.e0;
rfn=1.e0;
rfs=1.e0;
 
%     Set thermodynamic boundary conditions (for the seamount
%     problem, and other possible applications, lateral thermodynamic
%     boundary conditions are set equal to the initial conditions and
%     are held constant thereafter - users may, of course, create
%     variable boundary conditions):
for k=1:kbm1
	for j=1:jm
    	tbe(j,k)=tb(im,j,k);
        tbw(j,k)=tb(1,j,k);
        sbe(j,k)=sb(im,j,k);
        sbw(j,k)=sb(1,j,k);
    end
 
    for i=1:im
    	tbn(i,k)=tb(i,jm,k);
        tbs(i,k)=tb(i,1,k);
        sbn(i,k)=sb(i,jm,k);
        sbs(i,k)=sb(i,1,k);
    end
end

return 
end

