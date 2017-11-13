
function[dx,dy,cor,pdens,wetmask,...
    east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
    h,art,aru,arv,fsm,dum,dvm,...
    tsurf,ssurf,tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
    aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
    uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs] = wadseamount(dx,dy,cor,pdens,wetmask, ...
    east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
    h,hc,hco,hmin,hhi,nwad,time,art,aru,arv,fsm,dum,dvm,...
    tsurf,ssurf,tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
    aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
    uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,...
    e_atmos,aam,im,jm,kb,imm1,jmm1,kbm1,slmax,zz,tbias,sbias,grav,rhoref)

%  **********************************************************************
%  *                                                                    *
%  * FUNCTION    :  Sets up for seamount problem w/wad.                 *
%  *                                                                    *
%  **********************************************************************
%      real delh,delx,elejmid,elwjmid,ra,vel,hland  % lyo:%wad:define hland
%      integer nvar                                 % lyo:%wad:define nvar
%      integer i,j,k
%
%      Set delh > 1.0 for an island or delh < 1.0 for a seamount:
%
      delh=1.15e0;   % tne:%wad: island (for wad) =0.9e0 for orig.seamount
%
%      Grid size:
%
      delx=8000.e0;
%
% lyo:%wad:Define variable or uniform grid option:
      nvar=1;  % =1 for variable grid; =0 otherwise
% lyo:%wad:Special test case for smaller delx=4km:
      if ((nvar==0)||(im==131&&jm&&99))
      	delx=delx*0.5;
      end
%
%      Radius island or seamount:
%
      ra=50000.e0;   % tne:%wad:
%
%      Current velocity:
%
      vel=0.2e0;
%
%      Set up grid dimensions, areas of free surface cells, and
%      Coriolis parameter:
%
      for j=1:jm
        for i=1:im
%
%  For constant grid size:
%
%          dx(i,j)=delx
%          dy(i,j)=delx
%
%      For variable grid size:
%
%  lyo:%wad:variable or uniform grid--
          dx(i,j)=delx-(nvar)*delx*sin(pi*(i)/(im))/2.e0;
          dy(i,j)=delx-(nvar)*delx*sin(pi*(j)/(jm))/2.e0;
          cor(i,j)=1.e-4;
        end
      end
%{
%      Calculate horizontal coordinates of grid points and rotation
%      angle.
%
%      NOTE that this is introduced solely for the benefit of any post-
%      processing software, and in order to conform with the requirements
%      of the NetCDF Climate and Forecast (CF) Metadata Conventions.
%
%      There are four horizontal coordinate systems, denoted by the
%      subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
%      "e" is an elevation point and "c" is a cell corner), as shown
%      below. In addition, "east_*" is an easting and "north_*" is a
%      northing. Hence the coordinates of the "u" points are given by
%      (east_u,north_u).
%
%      Also, if the centre point of the cell shown below is at
%      (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
%      the coordinates of the western of the two "u" points,
%      (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
%      the two "v" points, and (east_c(i,j),north_c(i,j)) are the
%      coordinates of the southwestern corner point of the cell. The
%      southwest corner of the entire grid is at
%      (east_c(1,1),north_c(1,1)).
%
%                       |              |
%                     --c------v-------c--
%                       |              |
%                       |              |
%                       |              |
%                       |              |
%                       u      e       u
%                       |              |
%                       |              |
%                       |              |
%                       |              |
%                     --c------v-------c--
%                       |              |
%
%
%      NOTE that the following calculation of east_c and north_c only
%      works properly for a rectangular grid with east and north aligned
%      with i and j, respectively:
%}
      for j=1:jm
        east_c(1,j)=0.e0;
        for i=2:im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j);
        end
      end
%
      for i=1:im
        north_c(i,1)=0.e0;
        for j=2:jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1);
        end
      end
%
%      The following works properly for any grid:
%
%      Elevation points:
%
      for j=1:jm-1
        for i=1:im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)                          ...
                        +east_c(i,j+1)+east_c(i+1,j+1))/4.e0;
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)                       ...
                         +north_c(i,j+1)+north_c(i+1,j+1))/4.e0;
        end
      end
%
%      Extrapolate ends:
%
      for i=1:im-1
        east_e(i,jm)                                                      ...
          =((east_c(i,jm)+east_c(i+1,jm))*3.e0                            ...
             -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0;
        north_e(i,jm)                                                     ...
          =((north_c(i,jm)+north_c(i+1,jm))*3.e0                          ...
             -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0;
      end
%
      for j=1:jm-1
        east_e(im,j)                                                      ...
          =((east_c(im,j)+east_c(im,j+1))*3.e0                            ...
             -east_c(im-1,j)-east_c(im-1,j+1))/4.e0;
        north_e(im,j)                                                     ...
          =((north_c(im,j)+north_c(im,j+1))*3.e0                          ...
             -north_c(im-1,j)-north_c(im-1,j+1))/4.e0;
      end
%
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)                       ...
                     -(east_e(im-2:jm)+east_e(im,jm-2))/2.e0;
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)                    ...
                     -(north_e(im-2:jm)+north_e(im,jm-2))/2.e0;
%
%      u-points:
%
      for j=1:jm-1
        for i=1:im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0;
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0;
        end
      end
%
%      Extrapolate ends:
%
      for i=1:im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0;
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0;
      end
%
%      v-points:
%
      for j=1:jm
        for i=1:im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0;
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0;
        end
      end
%
%      Extrapolate ends:
%
      for j=1:jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0;
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0;
      end
%
%      rot is the angle (radians, anticlockwise) of the i-axis relative
%      to east, averaged to a cell centre:
%
%      (NOTE that the following calculation of rot only works properly
%      for this particular rectangular grid)
%
      for j=1:jm
        for i=1:im
          rot(i,j)=0.e0;
        end
      end
%
%      Define depth:
%
% tne:%wad:%lyo:%wad:
%      Note that unlike original seamount, h<0 is now water below MSL, and
%      h>0 is above the MSL and is either land (if h>hland, see below)
%      or is potential wet-and-dry region
%
      hmax=50.e0; hland=5.e0; %note all pnts are potential WAD if we set
                             %hland >= -hmax*(1.e0-delh) (=7.5m);
                             %i.e. "if" below is NOT satisfied
%
      for i=1:im
        for j=1:jm
%
          h(i,j)=-hmax*(1.e0-delh*exp(-((east_c(i,j)                       ...
                                         -east_c((im+1)/2,j))^2          ...
                                      +(north_c(i,j)                      ...
                                         -north_c(i,(jm+1)/2))^2)        ...
                                      /ra^2));
% lyo:%wad:Make region near top of seamount absolute land region:
          if(h(i,j)>hland)
          	h(i,j)=hhi+1.e0;
          end
%
        end
      end
%
%      Close the north and south boundaries to form a channel:
%
      for i=1:im
        h(i,1)=hhi+1.e0  ; % tne:%wad:
        h(i,jm)=hhi+1.e0 ; % tne:%wad:
      end
%
%      Calculate areas and masks:
%
% lyo:%wad:Define "h" consistent with WAD and calculate wetmask, cell areas
%          and fsm etc.
%      call areas_masks
      [art,aru,arv,fsm,dum,dvm,wetmask,h]=wadh(art,aru,arv,fsm,dum,dvm,wetmask, ...
        im,jm,dx,dy,h,hc,hco,hhi,hmin,nwad,slmax,time);
%
%      Adjust bottom topography so that cell to cell variations
%      in h for not exceed parameter slmax:
%
%      if(slmax<1.e0) call slpmax  %lyo:wad:now done in wadh above
%
%      Set initial conditions:
%
      for k=1:kbm1
        for j=1:jm
          for i=1:im
% lyo:wad:
%            tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
             if (h(i,j)>hhi)  %lyo:wad:stratified water below MSL:
                tb(i,j,k)=5.e0+15.e0*exp(zz(k)*(h(i,j)-hhi)/1000.e0);
             else                    %lyo:wad:well-mixed water above MSL:
                tb(i,j,k)=5.e0+15.e0;
             end
            tb(i,j,k)=tb(i,j,k)-tbias;
%
            sb(i,j,k)=35.e0-sbias ;
            tclim(i,j,k)=tb(i,j,k);
            sclim(i,j,k)=sb(i,j,k);
            ub(i,j,k)=vel*dum(i,j);
          end
        end
      end
%
%      Initialise uab and vab as necessary
%      (NOTE that these have already been initialised to zero in the
%      main program):
%
      for j=1:jm
        for i=1:im
          uab(i,j)=vel*dum(i,j);
        end
      end
%
%      Set surface boundary conditions, e_atmos, vflux, wusurf,
%      wvsurf, wtsurf, wssurf and swrad, as necessary
%      (NOTE:
%       1. These have all been initialised to zero in the main program.
%       2. The temperature and salinity of inflowing water must be
%          defined relative to tbias and sbias.):
%
      for j=1:jm
        for i=1:im
%      No conditions necessary for this problem
%
% lyo:%wad:%pom2k_bug:tsurf and ssurf were never defined, but should be:
             tsurf(i,j)=tb(i,j,1);
             ssurf(i,j)=sb(i,j,1);
%
        end
      end
%
%      Initialise elb, etb, dt and aam2d:
%
      for j=1:jm
        for i=1:im
          elb(i,j)=(-e_atmos(i,j)-hhi)*fsm(i,j);      % lyo:%wad:
% lyo:%wad:note:The following "if" is not satisfied if nwad=0 since     %
%      then fsm=wetmask at all cells                                    %
%      if (fsm(i,j)~=0.0&&wetmask(i,j)==0.0) % dry but not land
      		if (fsm(i,j)~=wetmask(i,j))                % dry but not land     ...
      			elb(i,j)=-h(i,j)+hco ; % note slightly smaller hco<hc to ensure that
      		  	                     % initially, cell is wet with elb+h=hco < hc
      		  etb(i,j)=elb(i,j)       ;        %lyo:%wad:
      		  dt(i,j)=h(i,j)+elb(i,j) ;        %lyo:%wad:
      		  aam2d(i,j)=aam(i,j,1)   ;
      		end
      	end
      end
%                                                                       %
% ----------------------------------------------------------------------%
% lyo:%wad:Make sure that initial wetmask remains compatible with       %
%      initial elb.  Also define marsh evaporation, +ve upward or       %
%      evaporate.  Also initialize wriv for rivers.                     %
% lyo:%wad:note:wriv is added here but the implementation following     %
%      Oey [1996; JPO, 26, 145-175; but see p.154 in particular] is not %
%      complete in this code.  E-mail me if want to implement it.       %
%                                                                       %
      if (nwad==1)
      	for j=1:jm;
      		for i=1:im
      			wetmask(i,j)=fsm(i,j);
      			if ((h(i,j)+elb(i,j))<=hc)
      				wetmask(i,j)=0.0;
      			end
      		end;
      	end
      end
%                                                                       %
      for j=1:jm;
      	for i=1:im
		      wmarsh(i,j)=0.0; wriv(i,j)=0.0;
		      if (fsm(i,j)~=wetmask(i,j))            % dry but not land
%      2 examples of wmarsh:
% --   wmarsh(i,j)=fsm(i,j)*10./art(i,j)
% --   wmarsh(i,j)=fsm(i,j)*(4.e-7)  %Gill's [1982] value
      		end
      	end;
      end
%                                                                       %
% ----------------------------------------------------------------------%
% lyo:%wad: Set up pdens before 1st call dens; used also in profq:      %
      for k=1:kbm1; for j=1:jm; for i=1:im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5;
      end; end; end
% ----------------------------------------------------------------------%
%                                                                       %
%      call dens(sb,tb,rho)
	[sb,tb,rho]=dens(sb,tb,rho,pdens, ...
	    im,jm,kbm1,tbias,sbias,grav,rhoref,zz,h,fsm);

%
%      Generated horizontally averaged density field (in this
%      application, the initial condition for density is a function
%      of z (the vertical cartesian coordinate) -- when this is not
%      so, make sure that rmean has been area averaged BEFORE transfer
%      to sigma coordinates):
%
      for k=1:kbm1
        for j=1:jm
          for i=1:im
            rmean(i,j,k)=rho(i,j,k);
          end
        end
      end
%
%      Set lateral boundary conditions, for use in subroutine bcond
%      (in the seamount problem, the east and west boundaries are open,
%      while the south and north boundaries are closed through the
%      specification of the masks fsm, dum and dvm):
%
      rfe=1.e0;      rfw=1.e0;      rfn=1.e0;      rfs=1.e0;
%
      for j=2:jmm1
        uabw(j)=uab(2,j)   ;
        uabe(j)=uab(imm1,j);
%
%      Set geostrophically conditioned elevations at the boundaries:
%
% lyo:%wad:Note keep (temporarily) all el* defined wrt MSL - ele, elw,
%          eln & els were initialized to be =0 in MAIN; then adjust
%          them later (below) with hhi:
%
        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1);
        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)         ;
      end
%
%      Adjust boundary elevations so that they are zero in the middle
%      of the channel:
%
      elejmid=ele(jmm1/2);
      elwjmid=elw(jmm1/2);
      for j=2:jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j);
        elw(j)=(elw(j)-elwjmid)*fsm(2,j) ;
      end
%
% tne:%wad:%lyo:wad:
      for i=1:im
         eln(i)=(0.-hhi)*fsm(i,jm);
         els(i)=(0.-hhi)*fsm(i,2) ;
      end
      
      for j=1:jm
         ele(j)=(ele(j)-hhi)*fsm(im,j);
         elw(j)=(elw(j)-hhi)*fsm(2,j) ;
      end
%
%      Set thermodynamic boundary conditions (for the seamount
%      problem, and other possible applications, lateral thermodynamic
%      boundary conditions are set equal to the initial conditions and
%      are held constant thereafter - users may, of course, create
%      variable boundary conditions):
%
      for k=1:kbm1
%
        for j=1:jm
          tbe(j,k)=tb(im,j,k);
          tbw(j,k)=tb(1,j,k) ;
          sbe(j,k)=sb(im,j,k);
          sbw(j,k)=sb(1,j,k) ;
        end
%
        for i=1:im
          tbn(i,k)=tb(i,jm,k);
          tbs(i,k)=tb(i,1,k) ;
          sbn(i,k)=sb(i,jm,k);
          sbs(i,k)=sb(i,1,k) ;
        end
%
      end
%

end
