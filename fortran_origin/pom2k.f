      program pom2k
C
C **********************************************************************
C *                                                                    *
C *   The last code change as rcorded in pom2k.change was on           *
C *                                                                    *
C *                     2006-05-03                                     *
C *                                  (adding IC from file)             *
C *                                                                    *
C * FUNCTION    :  This is a version of the three dimensional, time    *
C *                dependent, primitive equation, ocean model          *
C *                developed by Alan Blumberg and George Mellor with   *
C *                subsequent contributions by Leo Oey, Steve Brenner  *
C *                and others. It is now called the Princeton Ocean    *
C *                Model. Two references are:                          *
C *                                                                    *
C *                Blumberg, A.F. and G.L. Mellor; Diagnostic and      *
C *                  prognostic numerical circulation studies of the   *
C *                  South Atlantic Bight, J. Geophys. Res. 88,        *
C *                  4579-4592, 1983.                                  *
C *                                                                    *
C *                Blumberg, A.F. and G.L. Mellor; A description of a  *
C *                  three-dimensional coastal ocean circulation model,*
C *                  Three-Dimensional Coastal Ocean Models, Coastal   *
C *                  and Estuarine Sciences, 4, N.S. Heaps, ed.,       *
C *                  American Geophysical Union, 1-16, 1987.           *
C *                                                                    *
C *                In subroutine profq the model makes use of the      *
C *                turbulence closure sub-model described in:          *
C *                                                                    *
C *                Mellor, G.L. and T. Yamada; Development of a        *
C *                  turbulence closure model for geophysical fluid    *
C *                  problems, Rev. Geophys. Space Phys., 20, No. 4,   *
C *                  851-875, 1982.                                    *
C *            (note recent profq that includes breaking waves)        *
C *                                                                    *
C *                A user's guide is available:                        *
C *                                                                    *
C *                Mellor, G.L.; User's guide for a three-dimensional, *
C *                  primitive equation, numerical ocean model.        *
C *                  Princeton University Report, 1998.                *
C *                                                                    *
C *                In October 2001, the source code underwent          *
C *                revision by John Hunter of the University of        *
C *                Tasmania. Major aspects of the revision were:       *
C *                                                                    *
C *                (1) The revision was based on pom98 updated to      *
C *                    12/9/2001.                                      *
C *                (2) Declaration of all variables.                   *
C *                (3) Rationalisation of the input of all constants.  *
C *                (4) Modifications to the "printer" output.          *
C *                (5) Output to a netCDF file.                        *
C *                (6) Inclusion of surface freshwater flux.           *
C *                (7) Inclusion of atmospheric pressure.              *
C *                (8) Inclusion of an additional problem to check (6) *
C *                    and (7), above.                                 *
C *                (9) Inclusion of option for Smolarkiewicz           *
C *                    advection scheme.                               *
C *                                                                    *
C *                This revised version is functionally almost         *
C *                equivalent to pom98. The output to device 6 from    *
C *                the "seamount" problem should be almost the same,   *
C *                any differences being due to minor format changes   *
C *                and improvements in rounding.                       *
C *                                                                    *
C *                This revision was helped by the following people:   *
C *                Tal Ezer, Peter Holloway, George Mellor, Rich       *
C *                Signell, Ian Webster, Brian Williams and Emma Young.*
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *                                  GENERAL NOTES                     *
C *                                                                    *
C *                1. All units are S.I. (M.K.S.) unless otherwise     *
C *                   stated. NOTE that time is in days from the start *
C *                   of the run.                                      *
C *                                                                    *
C *                2. "b", <nothing> and "f" refers to backward,       *
C *                   central and forward time levels.                 *
C *                                                                    *
C *                3. NetCDF output may be used. In order to omit/use  *
C *                   netCDF, comment/uncomment all statements         *
C *                   carrying the comment "*netCDF*" at the end of    *
C *                   the line (or set netcdf_file='nonetcdf')         *
C *                                                                    *
C *                4. NetCDF is version 3. An attempt has been made to *
C *                   conform to the NetCDF Climate and Forecast (CF)  *
C *                   Metadata Conventions, but this may not yet be    *
C *                   complete (see:                                   *
C *                                                                    *
C *          http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html) *
C *                                                                    *
C *                5. In order to use netCDF, the program should be    *
C *                   compiled with the appropriate library. For       *
C *                   example, if using g77, you may need to type:     *
C *                                                                    *
C *                     g77 -o pom2k pom2k.f /usr/lib/libnetcdf.a      *
C *                                                                    *
C *                   You should also have the "include" file of       *
C *                   netCDF subroutines (pom2k.n).                    * 
C *                                                                    *
C *                6. In order to use netCDF, you may need to change   *
C *                   the name of the "include" file in the statement: *
C *                                                                    *
C *                     include '/usr/include/netcdf.inc'              *
C *                                                                    *
C *                   in subroutine write_netcdf                       *
C *                                                                    *
C **********************************************************************
C *                                                                    *
C *                                SOFTWARE LICENSING                  *
C *                                                                    *
C *                This program is free software; you can redistribute *
C *                it and/or modify it under the terms of the GNU      *
C *                General Public License as published by the Free     *
C *                Software Foundation, either Version 2 of the        *
C *                license, or (at your option) any later version.     *
C *                                                                    *
C *                This program is distributed in the hope that it     *
C *                will be useful, but without any warranty; without   *
C *                even the implied warranty of merchantability or     *
C *                fitness for a particular purpose. See the GNU       *
C *                General Public License for more details.            *
C *                                                                    *
C *                A copy of the GNU General Public License is         *
C *                available at http://www.gnu.org/copyleft/gpl.html   *
C *                or by writing to the Free Software Foundation, Inc.,*
C *                59 Temple Place - Suite 330, Boston, MA 02111, USA. *
C *                                                                    *
C **********************************************************************
C
C	  use netcdf
      implicit none
C
      include 'pom2k.c'
C
C     New declarations plus ispi,isp2i:
C
      real aam_init,atot
      real cbcmax,cbcmin,darea
      real days,dte2,dvol
      real eaver
      real horcon
      real ispi,isp2i
      real period,prtd1,prtd2
      real saver,smoth,sw,swtch
      real taver,tim0
      real vamax,vtot,tsalt
      real z0b
      real tatm,satm
      integer io(100),jo(100),ko(100)
      integer i,iend,iext,imax,ispadv,isplit,iswtch
      integer j,jmax
      integer k
      integer nadv,nbct,nbcs,nitera,nread
      integer iproblem
      logical lramp
      character*120 netcdf_file
C
C***********************************************************************
C
C     source should agree with source_c in pom2k.c and source_n in pom2k.n.
C
      source='pom2k  2006-05-03'
C
c     if(source.ne.source_c) then
c       write(6,7)
c   7   format(/'Incompatible versions of program and include files ',
c    $          '..... program terminated'/)
c       stop
c     endif
C
C***********************************************************************
C
      small=1.e-9           ! Small value
C
      pi=atan(1.0)*4.0    ! PI
C
C***********************************************************************
C
C     Input of filenames and constants:
C
C     NOTE that the array sizes im, jm and kb should be set in pom2k.c
C
C-----------------------------------------------------------------------
C
      title='Run 1                                   ' ! run's title
C
C-----------------------------------------------------------------------
C
      netcdf_file='pom2k.nc'  ! netCDF output file
c     netcdf_file='nonetcdf'  ! disable netCDF output
C
C-----------------------------------------------------------------------
C
C     Problem number:
C
C     iproblem      problem      initialisation
C                    type          subroutine
C
C         1        seamount       seamount
C
C         2        conservation   box
C                  box
C
C         3        IC from file   file2ic
C
      iproblem=1
C
C-----------------------------------------------------------------------
C
C       mode                     description
C
C        2        2-D calculation (bottom stress calculated in advave)
C
C        3        3-D calculation (bottom stress calculated in profu,v)
C
C        4        3-D calculation with t and s held fixed
C
      mode=3
C
C-----------------------------------------------------------------------
C
C     Advection scheme:
C
C      nadv     Advection scheme
C
C        1       Centred scheme, as originally provide in POM
C        2       Smolarkiewicz iterative upstream scheme, based on
C                subroutines provided by Gianmaria Sannino and Vincenzo
C                Artale
C
      nadv=1
C
C-----------------------------------------------------------------------
C
C     Constants for Smolarkiewicz iterative upstream scheme.
C
C     Number of iterations. This should be in the range 1 - 4. 1 is
C     standard upstream differencing; 3 adds 50% CPU time to POM:
C
      nitera=2
C
C     Smoothing parameter. This should preferably be 1, but 0 < sw < 1
C     gives smoother solutions with less overshoot when nitera > 1:
C
      sw=0.50
C
C-----------------------------------------------------------------------
C
C     Index to indicate whether run to start from restart file
C     (nread=0: no restart input file; nread=1: restart input file):
C
      nread=0
C
C-----------------------------------------------------------------------
C
C     External (2-D) time step (secs.) according to CFL:
C
      dte=6.0
C
C-----------------------------------------------------------------------
C
C     <Internal (3-D) time step>/<External (2-D) time step>
C     (dti/dte; dimensionless):
C
      isplit=30
C
C-----------------------------------------------------------------------
C
C     Date and time of start of initial run of model in format (i.e.
C     UDUNITS convention)
C
C       YYYY-MM-DD HH:MM:SS <+/->HH:MM
C
C     where "<+/->HH:MM" is the time zone (positive eastwards from
C     Coordinated Universal Time). NOTE that the climatological time
C     axis (i.e. beginning of year zero, which does not exist in the
C     real-world calendar) has been used here. Insert your own date
C     and time as required:
C
      time_start='2000-01-01 00:00:00 +00:00'
C
C-----------------------------------------------------------------------
C
      days=0.025       ! run duration in days
C
C-----------------------------------------------------------------------
C
      prtd1=0.0125     ! Initial print interval (days)
C
C-----------------------------------------------------------------------
C
      prtd2=1.0         ! Final print interval (days)
C
C-----------------------------------------------------------------------
C
      swtch=1000.0      ! Time to switch from prtd1 to prtd2 
C
C-----------------------------------------------------------------------
C
      iskp=4             ! Printout skip interval in i 
C
C-----------------------------------------------------------------------
C
      jskp=3             ! Printout skip interval in j
C
C-----------------------------------------------------------------------
C
C     Logical for inertial ramp (.true. if inertial ramp to be applied
C     to wind stress and baroclinic forcing, otherwise .false.)
C
      lramp=.false.
C
C-----------------------------------------------------------------------
C
C     Reference density (recommended values: 1025 for seawater,
C     1000 for freswater; S.I. units):
C
      rhoref=1025.0
C
C-----------------------------------------------------------------------
C
      tbias=0.0         ! Temperature bias (deg. C)
C
C-----------------------------------------------------------------------
C
      sbias=0.0         ! Salinity bias
C
C-----------------------------------------------------------------------
C
      grav=9.8060       ! gravity constant (S.I. units)
C
C-----------------------------------------------------------------------
C
      kappa=0.40        ! von Karman's constant
C
C-----------------------------------------------------------------------
C
      z0b=.010          ! Bottom roughness (metres)
C
C-----------------------------------------------------------------------
C
      cbcmin=.00250     ! Minimum bottom friction coeff.
C
C-----------------------------------------------------------------------
C
      cbcmax=1.0        ! Maximum bottom friction coeff.
C
C-----------------------------------------------------------------------
C
      horcon=0.20       ! Smagorinsky diffusivity coeff.
C
C-----------------------------------------------------------------------
C
C     Inverse horizontal turbulent Prandtl number
C     (ah/am; dimensionless):
C
C     NOTE that tprni=0.0 yields zero horizontal diffusivity!
C
      tprni=.20
C
C-----------------------------------------------------------------------
C
C     Background viscosity used in subroutines profq, proft, profu and
C     profv (S.I. units):
C
      umol=2.e-5
C
C-----------------------------------------------------------------------
C
C     Maximum depth used in radiation boundary condition in subroutine
C     bcond (metres):
C
      hmax=4500.0
C
C-----------------------------------------------------------------------
C
C     Maximum magnitude of vaf (used in check that essentially tests
C     for CFL violation):
C
      vmaxl=100.0
C
C-----------------------------------------------------------------------
C
C     Maximum allowable value of:
C
C       <difference of depths>/<sum of depths>
C
C     for two adjacent cells (dimensionless). This is used in subroutine
C     slpmax. If >= 1, then slpmax is not applied:
C
      slmax=2.0
C
C-----------------------------------------------------------------------
C
C     Integers defining the number of logarithmic layers at the   
C     surface and bottom (used by subroutine depth). The number of
C     logarithmic layers are kl1-2 at the surface and kb-kl2-1
C     at the bottom. For no log portions, set kl1=2 and kl2=kb-1:
C
      kl1=6
      kl2=kb-2
C
C-----------------------------------------------------------------------
C
C     Water type, used in subroutine proft.
C
C       ntp    Jerlov water type
C
C        1            i
C        2            ia
C        3            ib
C        4            ii
C        5            iii
C
      ntp=2
C
C-----------------------------------------------------------------------
C
C     Surface temperature boundary condition, used in subroutine proft:
C
C       nbct   prescribed    prescribed   short wave
C              temperature      flux      penetration
C
C        1        no           yes           no
C        2        no           yes           yes
C        3        yes          no            no
C        4        yes          no            yes
C
      nbct=1
C
C-----------------------------------------------------------------------
C
C     Surface salinity boundary condition, used in subroutine proft:
C
C       nbcs   prescribed    prescribed
C               salinity      flux
C
C        1        no           yes
C        3        yes          no
C
C     NOTE that only 1 and 3 are allowed for salinity.
C
      nbcs=1
C
C-----------------------------------------------------------------------
C
C     Step interval during which external (2-D) mode advective terms are
C     not updated (dimensionless):
C
      ispadv=5
C
C-----------------------------------------------------------------------
C
C     Constant in temporal filter used to prevent solution splitting
C     (dimensionless):
C
      smoth=0.100
C
C-----------------------------------------------------------------------
C
C     Weight used for surface slope term in external (2-D) dynamic
C     equation (a value of alpha = 0.0 is perfectly acceptable, but the
C     value, alpha=.2250 permits a longer time step):
C
      alpha=0.2250
C
C-----------------------------------------------------------------------
C
C     Initial value of aam:
C
      aam_init=500.0
C
C     End of input of constants
C***********************************************************************
C
C --- Above are the default parameters, alternatively one can 
C --- use parameters from a file created by runscript runpom2k
C
      include 'params'
C
C***********************************************************************
C
C     Calculate some constants:
C
      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2
C
      iend=max0(nint(days*24.0*3600.0/dti),2)
      iprint=nint(prtd1*24.0*3600.0/dti)
      iswtch=nint(swtch*24.0*3600.0/dti)
C
      ispi=1.0/float(isplit)
      isp2i=1.0/(2.0*float(isplit))
C
C-----------------------------------------------------------------------
C
C     Print initial summary:
C
      write(6,'(/,'' source   = '',a40)') source
      write(6,'('' title      = '',a40/)') title
      write(6,'('' iproblem   = '',i10)') iproblem 
      write(6,'('' mode       = '',i10)') mode
      write(6,'('' nadv       = '',i10)') nadv
      write(6,'('' nitera     = '',i10)') nitera
      write(6,'('' sw         = '',f10.4)') sw
      write(6,'('' nread      = '',i10)') nread   
      write(6,'('' dte        = '',f10.2)') dte
      write(6,'('' dti        = '',f10.1)') dti
      write(6,'('' isplit     = '',i10)') isplit 
      write(6,'('' time_start = '',a26)') time_start
      write(6,'('' days       = '',f10.4)') days
      write(6,'('' iend       = '',i10)') iend
      write(6,'('' prtd1      = '',f10.4)') prtd1
      write(6,'('' iprint     = '',i10)') iprint 
      write(6,'('' prtd2      = '',f10.4)') prtd2
      write(6,'('' swtch      = '',f10.2)') swtch 
      write(6,'('' iswtch     = '',i10)') iswtch 
      write(6,'('' iskp, jskp = '',i5'','',i5)') iskp,jskp   
      write(6,'('' lramp      = '',l10)') lramp  
      write(6,'('' rhoref     = '',f10.3)') rhoref
      write(6,'('' tbias      = '',f10.3)') tbias 
      write(6,'('' sbias      = '',f10.3)') sbias 
      write(6,'('' grav       = '',f10.4)') grav  
      write(6,'('' kappa      = '',f10.4)') kappa  
      write(6,'('' z0b        = '',f10.6)') z0b         
      write(6,'('' cbcmin     = '',f10.6)') cbcmin      
      write(6,'('' cbcmax     = '',f10.6)') cbcmax      
      write(6,'('' horcon     = '',f10.3)') horcon      
      write(6,'('' tprni      = '',f10.4)') tprni      
      write(6,'('' umol       = '',f10.4)') umol       
      write(6,'('' hmax       = '',f10.2)') hmax       
      write(6,'('' vmaxl      = '',f10.4)') vmaxl      
      write(6,'('' slmax      = '',f10.4)') slmax      
      write(6,'('' kl1, kl2   = '',i5,'','',i5)') kl1,kl2
      write(6,'('' ntp        = '',i10)') ntp  
      write(6,'('' nbct       = '',i10)') nbct  
      write(6,'('' nbcs       = '',i10)') nbcs  
      write(6,'('' ispadv     = '',i10)') ispadv
      write(6,'('' smoth      = '',f10.4)') smoth     
      write(6,'('' alpha      = '',f10.4)') alpha     
C
C-----------------------------------------------------------------------
C
C     Initialise boundary arrays:
C
      do i=1,im
        vabn(i)=0.0
        vabs(i)=0.0
        eln(i)=0.0
        els(i)=0.0
        do k=1,kb
          vbn(i,k)=0.0
          vbs(i,k)=0.0
          tbn(i,k)=0.0
          tbs(i,k)=0.0
          sbn(i,k)=0.0
          sbs(i,k)=0.0
        end do
      end do
C
      do j=1,jm
        uabe(j)=0.0
        uabw(j)=0.0
        ele(j)=0.0
        elw(j)=0.0
        do k=1,kb
          ube(j,k)=0.0
          ubw(j,k)=0.0
          tbe(j,k)=0.0
          tbw(j,k)=0.0
          sbe(j,k)=0.0
          sbw(j,k)=0.0
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     Initialise 2-D and 3-D arrays for safety (this may be overwritten
C     later):
C
      do j=1,jm
        do i=1,im
          uab(i,j)=0.0
          vab(i,j)=0.0
          elb(i,j)=0.0
          etb(i,j)=0.0
          e_atmos(i,j)=0.0
          vfluxb(i,j)=0.0
          vfluxf(i,j)=0.0
          wusurf(i,j)=0.0
          wvsurf(i,j)=0.0
          wtsurf(i,j)=0.0
          wssurf(i,j)=0.0
          swrad(i,j)=0.0
          drx2d(i,j)=0.0
          dry2d(i,j)=0.0
        end do
      end do
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            ub(i,j,k)=0.0
            vb(i,j,k)=0.0
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     Set up sigma layers:
C
      if(iproblem.ne.3) call depth
C
C-----------------------------------------------------------------------
C
C     Read in grid data, and initial and lateral boundary conditions:
C
      if(iproblem.eq.1) then
        call seamount
      else if(iproblem.eq.2) then
        call box
      else if(iproblem.eq.3) then
        call file2ic
      else
        write(6,8)
    8   format(/' Invalid value of iproblem ..... program terminated'/)
        stop
      endif
C
C     Inertial period for temporal filter:
C
      period=(2.0*pi)/abs(cor(im/2,jm/2))/86400.0
C
C     Initialise time:
C
      tim0=0.0
      time=0.0
C
C     Initial conditions:
C
C     NOTE that lateral thermodynamic boundary conditions are often set
C     equal to the initial conditions and are held constant thereafter.
C     Users can of course create variable boundary conditions.
C
      do i=1,im
        do j=1,jm
          ua(i,j)=uab(i,j)
          va(i,j)=vab(i,j)
          el(i,j)=elb(i,j)
          et(i,j)=etb(i,j)
          etf(i,j)=et(i,j)
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          w(i,j,1)=vfluxf(i,j)
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            l(i,j,k)=0.1*dt(i,j)
            q2b(i,j,k)=small
            q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k)
            kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k))
            km(i,j,k)=kh(i,j,k)
            kq(i,j,k)=kh(i,j,k)
            aam(i,j,k)=aam_init
          end do
        end do
      end do
C
      do k=1,kbm1
        do i=1,im
          do j=1,jm
            q2(i,j,k)=q2b(i,j,k)
            q2l(i,j,k)=q2lb(i,j,k)
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            u(i,j,k)=ub(i,j,k)
            v(i,j,k)=vb(i,j,k)
          end do
        end do
      end do
C
      call dens(s,t,rho)
C
      call baropg
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do
C
C     Calculate bottom friction coefficient:
C
      do j=1,jm
        do i=1,im
          cbc(i,j)=(kappa/log((1.0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
C
C     If the following is invoked, then it is probable that the wrong
C     choice of z0b or vertical spacing has been made:
C
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
C
C     Calculate external (2-D) CFL time step:
C
      do j=1,jm
        do i=1,im
          tps(i,j)=0.50/sqrt(1.0/dx(i,j)**2+1.0/dy(i,j)**2)
     $               /sqrt(grav*(h(i,j)+small))*fsm(i,j)
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following data are needed for a seamless restart. if nread=1,
C     data had been created by a previous run (see write(71) at end of
C     this program). nread=0 denotes a first time run.
C
      if(nread.eq.1)
     $  read(70) tim0,
     $           wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $           utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,
     $           adx2d,ady2d,advua,advva,
     $           km,kh,kq,l,q2,q2b,aam,q2l,q2lb
C
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do
C
      time=tim0
C
C-----------------------------------------------------------------------
C
C     Print geometry and other initial fields (select statements as
C     desired):
C
      call prxy('grid increment in x, dx                 ',
     $          time,dx ,im,iskp,jm,jskp,0.0)
C
      call prxy('grid increment in y, dy                 ',
     $          time,dy ,im,iskp,jm,jskp,0.0)
C
      call prxy('Easting of elevation points, east_e     ',
     $          time,east_e ,im,iskp,jm,jskp,0.0)
C
      call prxy('Northing of elevation points, north_e   ',
     $          time,north_e,im,iskp,jm,jskp,0.0)
C
      call prxy('Easting of cell corners, east_c         ',
     $          time,east_c ,im,iskp,jm,jskp,0.0)
C
      call prxy('Northing of cell corners, north_c       ',
     $          time,north_c,im,iskp,jm,jskp,0.0)
C
      call prxy('Rotation angle of x-axis wrt. east, rot ',
     $          time,rot,im,iskp,jm,jskp,0.0)
C
      call prxy('Undisturbed water depth, h              ',
     $          time,h  ,im,iskp,jm,jskp,1.e1)
C
      call prxy('Free surface mask, fsm                  ',
     $          time,fsm,im,iskp,jm,jskp,1.0)
C
      call prxy('u-velocity mask, dum                    ',
     $          time,dum,im,iskp,jm,jskp,1.0)
C
      call prxy('v-velocity mask, dvm                    ',
     $          time,dvm,im,iskp,jm,jskp,1.0)
C
      call prxy('External (2-D) CFL time step, tps       ',
     $          time,tps,im,iskp,jm,jskp,1.0)
C
C     Set sections for output:
C
      ko(1)=1
      ko(2)=kb/2
      ko(3)=kb-1
C
      call prxyz('Horizontally-averaged rho, rmean        ',
     $           time,rmean,im,iskp,jm,jskp,kb,ko,3,1.e-5)
C
C     Set sections for output:
C
      jo(1)=1
      jo(2)=jm/2
      jo(3)=jm-1
C
      call prxz('Horizontally-averaged rho, rmean        ',
     $          time,rmean,im,iskp,jm,kb,jo,3,1.e-5,dt,zz)
C
C     Set sections for output:
C
      io(1)=1
      io(2)=im/2
      io(3)=im-1
C
      call pryz('Horizontally-averaged rho, rmean        ',
     $          time,rmean,im,jm,jskp,kb,io,3,1.e-5,dt,zz)
C
C-----------------------------------------------------------------------
C
C     Initial conditions:
C
C     Select print statements in printall as desired:
C
      call printall
C
C-----------------------------------------------------------------------
C
C     Initialise netCDF output and output initial set of data:
C
        if(netcdf_file.ne.'nonetcdf') then
C      call write_netcdf(netcdf_file,1)                        ! *netCDF*
C      call write_netcdf(netcdf_file,2)                        ! *netCDF*
        endif
C
C-----------------------------------------------------------------------
C
      do 9000 iint=1,iend      !  Begin internal (3-D) mode
C
        time=dti*float(iint)/86400.0+tim0
C
        if(lramp) then
          ramp=time/period
          if(ramp.gt.1.0) ramp=1.0
        else
          ramp=1.0
        endif
C
C       write(6,2) mode,iint,time
C   2   format(' mode,iint,time =',2i5,f9.2)
C
C-----------------------------------------------------------------------
C
C     Set time dependent, surface and lateral boundary conditions.
C     The latter will be used in subroutine bcond. Users may
C     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
C     wssurf, swrad and vflux.
C
C     Introduce simple wind stress. Value is negative for westerly or
C     southerly winds. The following wind stress has been tapered
C     along the boundary to suppress numerically induced oscilations
C     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
C     To make a healthy surface Ekman layer, it would be well to set
C     kl1=9.
C
        do j=2,jmm1
          do i=2,imm1
c
      if(iproblem.ne.3) then     ! constant wind read in file2ic
c
c           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1))
            wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))
     $                    *.250*(dvm(i,j+1)+dvm(i-1,j+1)
     $                          +dvm(i-1,j)+dvm(i,j))
C --- no wind ----
c           wusurf(i,j)=0.0
            wvsurf(i,j)=0.0
       endif
            e_atmos(i,j)=0.0
            vfluxf(i,j)=0.0
C
C     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
C     the sea surface. See calculation of elf(i,j) below and subroutines
C     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
C     is no net flow across lateral boundaries, the basin volume will be
C     constant; if also vflux(i,j).ne.0, then, for example, the average
C     salinity will change and, unrealistically, so will total salt. 
C
            w(i,j,1)=vfluxf(i,j)
C
C     Set wtsurf to the sensible heat, the latent heat (which involves
C     only the evaporative component of vflux) and the long wave
C     radiation:
C
            wtsurf(i,j)=0.0
C
C     Set swrad to the short wave radiation:
C
            swrad(i,j)=0.0
C
C     To account for change in temperature of flow crossing the sea
C     surface (generally quite small compared to latent heat effect)
C
            tatm=t(i,j,1)+tbias    ! an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias)
C
C     Set the salinity of water vapor/precipitation which enters/leaves
C     the atmosphere (or e.g., an ice cover)
C    
            satm=0.0              
            wssurf(i,j)=            vfluxf(i,j)*(satm-s(i,j,1)-sbias)  
C
          end do
        end do
C
C-----------------------------------------------------------------------
C
C     Set lateral viscosity:
C
C     If mode=2 then initial values of aam2d are used. If one wishes
C     to use Smagorinsky lateral viscosity and diffusion for an
C     external (2-D) mode calculation, then appropiate code can be
C     adapted from that below and installed just before the end of the
C     "if(mode.eq.2)" loop in subroutine advave.
C
C     Calculate Smagorinsky lateral viscosity:
C
C       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
C                                     +.5*(du/dy+dv/dx)**2) )
C
        if(mode.ne.2) then
          call advct(a,c,ee)
          call baropg
C
          do k=1,kbm1
            do j=2,jmm1
              do i=2,imm1
                aam(i,j,k)=horcon*dx(i,j)*dy(i,j)
     $                      *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2
     $                            +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2
     $                      +.50*(.250*(u(i,j+1,k)+u(i+1,j+1,k)
     $                                   -u(i,j-1,k)-u(i+1,j-1,k))
     $                      /dy(i,j)
     $                      +.250*(v(i+1,j,k)+v(i+1,j+1,k)
     $                             -v(i-1,j,k)-v(i-1,j+1,k))
     $                      /dx(i,j)) **2)
              end do
            end do
          end do
C
C     Form vertical averages of 3-D fields for use in external (2-D)
C     mode:
C
          do j=1,jm
            do i=1,im
              adx2d(i,j)=0.0
              ady2d(i,j)=0.0
              drx2d(i,j)=0.0
              dry2d(i,j)=0.0
              aam2d(i,j)=0.0
            end do
          end do
C
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k)
                ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k)
                drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
                dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
                aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k)
              end do
            end do
          end do
C
          call advave(tps)
C
          do j=1,jm
            do i=1,im
              adx2d(i,j)=adx2d(i,j)-advua(i,j)
              ady2d(i,j)=ady2d(i,j)-advva(i,j)
            end do
          end do
C
        endif
C
        do j=1,jm
          do i=1,im
            egf(i,j)=el(i,j)*ispi
          end do
        end do
C
        do j=1,jm
          do i=2,im
            utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
          end do
        end do
        do j=2,jm
          do i=1,im
            vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
          end do
        end do
C
C-----------------------------------------------------------------------
C
        do 8000 iext=1,isplit    ! Begin external (2-D) mode 
C
C         write(6,3) iext,time
C   3     format(' iext,time =',i5,f9.2)
C
          do j=2,jm
            do i=2,im
              fluxua(i,j)=.250*(d(i,j)+d(i-1,j))
     $                     *(dy(i,j)+dy(i-1,j))*ua(i,j)
              fluxva(i,j)=.250*(d(i,j)+d(i,j-1))
     $                     *(dx(i,j)+dx(i,j-1))*va(i,j)
            end do
          end do
C
C     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
C     with pom98.f. See also modifications to subroutine vertvl.
C
          do j=2,jmm1
            do i=2,imm1
              elf(i,j)=elb(i,j)
     $                  +dte2*(-(fluxua(i+1,j)-fluxua(i,j)
     $                          +fluxva(i,j+1)-fluxva(i,j))/art(i,j)
     $                          -vfluxf(i,j))
            end do
          end do
C
          call bcond(1)

          if(mod(iext,ispadv).eq.0) call advave(tps)
C
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=adx2d(i,j)+advua(i,j)
     $                  -aru(i,j)*.250
     $                    *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))
     $                     +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
     $                  +.250*grav*(dy(i,j)+dy(i-1,j))
     $                    *(d(i,j)+d(i-1,j))
     $                    *((1.0-2.0*alpha)
     $                       *(el(i,j)-el(i-1,j))
     $                      +alpha*(elb(i,j)-elb(i-1,j)
     $                             +elf(i,j)-elf(i-1,j))
     $                      +e_atmos(i,j)-e_atmos(i-1,j))
     $                  +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
            end do
          end do
C
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))
     $                    *aru(i,j)*uab(i,j)
     $                  -4.0*dte*uaf(i,j))
     $                 /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))
     $                     *aru(i,j))
            end do
          end do
C
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=ady2d(i,j)+advva(i,j)
     $                  +arv(i,j)*.250
     $                    *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))
     $                     +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
     $                  +.250*grav*(dx(i,j)+dx(i,j-1))
     $                    *(d(i,j)+d(i,j-1))
     $                    *((1.0-2.0*alpha)*(el(i,j)-el(i,j-1))
     $                      +alpha*(elb(i,j)-elb(i,j-1)
     $                             +elf(i,j)-elf(i,j-1))
     $                      +e_atmos(i,j)-e_atmos(i,j-1))
     $                  +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
            end do
          end do
C
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))
     $                    *vab(i,j)*arv(i,j)
     $                  -4.0*dte*vaf(i,j))
     $                 /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))
     $                     *arv(i,j))
            end do
          end do
C
          call bcond(2)
C
          if(iext.eq.(isplit-2))then
            do j=1,jm
              do i=1,im
                etf(i,j)=.250*smoth*elf(i,j)
              end do
            end do
C
          else if(iext.eq.(isplit-1)) then
C
            do j=1,jm
              do i=1,im
                etf(i,j)=etf(i,j)+.50*(1.-.50*smoth)*elf(i,j)
              end do
            end do
C
          else if(iext.eq.isplit) then
C
            do j=1,jm
              do i=1,im
                etf(i,j)=(etf(i,j)+.50*elf(i,j))*fsm(i,j)
              end do
            end do
C
          endif
C
C     Stop if velocity condition violated (generally due to CFL
C     criterion not being satisfied):
C
          vamax=0.0
C
          do j=1,jm
            do i=1,im
              if(abs(vaf(i,j)).ge.vamax) then
                vamax=abs(vaf(i,j))
	        imax=i
	        jmax=j
              endif
            end do
          end do
C
          if(vamax.le.vmaxl) then
C
C     Apply filter to remove time split and reset time sequence:
C
            do j=1,jm
              do i=1,im
                ua(i,j)=ua(i,j)
     $                   +.50*smoth*(uab(i,j)-2.0*ua(i,j)+uaf(i,j))
                va(i,j)=va(i,j)
     $                   +.50*smoth*(vab(i,j)-2.0*va(i,j)+vaf(i,j))
                el(i,j)=el(i,j)
     $                   +.50*smoth*(elb(i,j)-2.0*el(i,j)+elf(i,j))
                elb(i,j)=el(i,j)
                el(i,j)=elf(i,j)
                d(i,j)=h(i,j)+el(i,j)
                uab(i,j)=ua(i,j)
                ua(i,j)=uaf(i,j)
                vab(i,j)=va(i,j)
                va(i,j)=vaf(i,j)
              end do
            end do
C
            if(iext.ne.isplit) then
              do j=1,jm
                do i=1,im
                  egf(i,j)=egf(i,j)+el(i,j)*ispi
                end do
              end do
              do j=1,jm
                do i=2,im
                  utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
                end do
              end do
              do j=2,jm
                do i=1,im
                  vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
                end do
              end do
            endif
C
          endif
C
 8000 continue        ! End of external (2-D) mode
C
C-----------------------------------------------------------------------
C
        if(vamax.le.vmaxl) then
C
C     Continue with internal (3-D) mode calculation:
C
          if((iint.ne.1.or.tim0.ne.0.0).and.mode.ne.2) then
C
C     Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=2,im
                  u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                     (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
                end do
              end do
            end do
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=2,jm
                do i=1,im
                  v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                     (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
                end do
              end do
            end do
C
C     vertvl calculates w from u, v, dt (h+et), etf and etb:
C
            call vertvl(a,c)
            call bcond(5)
C
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  uf(i,j,k)=0.0
                  vf(i,j,k)=0.0
                end do
              end do
            end do
C
C     Calculate q2f and q2lf using uf, vf, a and c as temporary
C     variables:
C
            call advq(q2b,q2,uf,a,c)
            call advq(q2lb,q2l,vf,a,c)
            call profq(a,c,tps,dtef)
            call bcond(6)
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  q2(i,j,k)=q2(i,j,k)
     $                       +.50*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                                    -2.0*q2(i,j,k))
                  q2l(i,j,k)=q2l(i,j,k)
     $                       +.50*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                                    -2.0*q2l(i,j,k))
                  q2b(i,j,k)=q2(i,j,k)
                  q2(i,j,k)=uf(i,j,k)
                  q2lb(i,j,k)=q2l(i,j,k)
                  q2l(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
C
C     Calculate tf and sf using uf, vf, a and c as temporary variables:
C
            if(mode.ne.4) then
C
              if(nadv.eq.1) then
C
                call advt1(tb,t,tclim,uf,a,c)
                call advt1(sb,s,sclim,vf,a,c)
C
              else if(nadv.eq.2) then
C
                call advt2(tb,t,tclim,uf,a,c,nitera,sw)
                call advt2(sb,s,sclim,vf,a,c,nitera,sw)
C
              else
C
                write(6,9)
    9           format(/'Invalid value for nadv ..... ',
     $                 'program terminated'/)
                stop
C
              endif
C
              call proft(uf,wtsurf,tsurf,nbct,tps)
              call proft(vf,wssurf,ssurf,nbcs,tps)
              call bcond(4)
C
              do k=1,kb
                do j=1,jm
                  do i=1,im
                    t(i,j,k)=t(i,j,k)
     $                        +.50*smoth*(uf(i,j,k)+tb(i,j,k)
     $                                     -2.0*t(i,j,k))
                    s(i,j,k)=s(i,j,k)
     $                        +.50*smoth*(vf(i,j,k)+sb(i,j,k)
     $                                     -2.0*s(i,j,k))
                    tb(i,j,k)=t(i,j,k)
                    t(i,j,k)=uf(i,j,k)
                    sb(i,j,k)=s(i,j,k)
                    s(i,j,k)=vf(i,j,k)
                  end do
                end do
              end do
C
              call dens(s,t,rho)
C
            endif
C
C     Calculate uf and vf:
C
            call advu
            call advv
            call profu
            call profv
            call bcond(3)
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)
     $                      +(uf(i,j,k)+ub(i,j,k)-2.0*u(i,j,k))*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  u(i,j,k)=u(i,j,k)
     $                      +.50*smoth*(uf(i,j,k)+ub(i,j,k)
     $                                   -2.0*u(i,j,k)-tps(i,j))
                end do
              end do
            end do
C
            do j=1,jm
              do i=1,im
                tps(i,j)=0.0
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)
     $                      +(vf(i,j,k)+vb(i,j,k)-2.0*v(i,j,k))*dz(k)
                end do
              end do
            end do
C
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  v(i,j,k)=v(i,j,k)
     $                      +.50*smoth*(vf(i,j,k)+vb(i,j,k)
     $                                   -2.0*v(i,j,k)-tps(i,j))
                end do
              end do
            end do
C
            do k=1,kb
              do j=1,jm
                do i=1,im
                  ub(i,j,k)=u(i,j,k)
                  u(i,j,k)=uf(i,j,k)
                  vb(i,j,k)=v(i,j,k)
                  v(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
C
          endif
C
          do j=1,jm
            do i=1,im
              egb(i,j)=egf(i,j)
              etb(i,j)=et(i,j)
              et(i,j)=etf(i,j)
              dt(i,j)=h(i,j)+et(i,j)
              utb(i,j)=utf(i,j)
              vtb(i,j)=vtf(i,j)
              vfluxb(i,j)=vfluxf(i,j)
            end do
          end do
C
        endif
C
C-----------------------------------------------------------------------
C
C     Beginning of print section:
C
        if(iint.ge.iswtch) iprint=nint(prtd2*24.0*3600.0/dti)
C
        if(mod(iint,iprint).eq.0.or.vamax.gt.vmaxl) then
C
          write(6,4) time,iint,iext,iprint
    4     format(/
     $    '**************************************************',
     $    '**************************************************',
     $    '*************************'//
     $    ' time =',f9.4,', iint =',i8,', iext =',i8,', iprint =',i8,//)
C
C     Select print statements in printall as desired:
C
          call printall
C
          vtot=0.0
          atot=0.0
          taver=0.0
          saver=0.0
          eaver=0.0
          do k=1,kbm1
            do j=1,jm
              do i=1,im
                darea=dx(i,j)*dy(i,j)*fsm(i,j)
                dvol=darea*dt(i,j)*dz(k)
                vtot=vtot+dvol
                taver=taver+tb(i,j,k)*dvol
                saver=saver+sb(i,j,k)*dvol
              end do
            end do
          end do
C
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              atot=atot+darea
              eaver=eaver+et(i,j)*darea
            end do
          end do
C
          taver=taver/vtot
          saver=saver/vtot
          eaver=eaver/atot
          tsalt=(saver+sbias)*vtot
C
          write(6,5) vtot,atot,eaver,taver,saver,tsalt
    5     format('vtot = ',e16.7,'   atot = ',e16.7,
     $           '  eaver =',e16.7/'taver =',e16.7,
     $           '   saver =',e16.7,'  tsalt =',e16.7)  
C
C     Write netCDF output:
C
            if(netcdf_file.ne.'nonetcdf') then
C          call write_netcdf(netcdf_file,2)                    ! *netCDF*
            endif
C
          if(vamax.gt.vmaxl) then
C
            write(6,4) time,iint,iext,iprint
C
            call printall
C
            write(6,6) vamax,imax,jmax
    6       format(///////////////////
     $             '************************************************'/
     $             '************ abnormal job end ******************'/
     $             '************* user terminated ******************'/
     $             '************************************************'/
     $             ' vamax =',e12.3,'   imax,jmax =',2i5)
C
C     Close netCDF file:
C
              if(netcdf_file.ne.'nonetcdf') then
C            call write_netcdf(netcdf_file,3)                  ! *netCDF*
              endif
C
            stop
C
          endif
C
        endif
C
C     End of print section
C
C-----------------------------------------------------------------------
C
 9000 continue       !  End of internal (3-D) mode
C
C-----------------------------------------------------------------------
C
      write(6,4) time,iint,iext,iprint
C
C     Set levels for output:
C
      ko(1)=1
      ko(2)=2
      ko(3)=kb/2
      ko(4)=kb-1
      ko(5)=kb
C
C     call prxyz('Vertical velocity, w                    ',
C    $           time,w       ,im,iskp,jm,jskp,kb,ko,5,-1.0)
C
C     call prxyz('Turbulent kinetic energy x 2, q2        ',
C    $           time,q2      ,im,iskp,jm,jskp,kb,ko,5,-1.0)
C
C     Save this data for a seamless restart:
C
      write(71) time,
     $  wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,
     $  utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,
     $  km,kh,kq,l,q2,q2b,aam,q2l,q2lb
C
C     Close netCDF file:
C
        if(netcdf_file.ne.'nonetcdf') then
C      call write_netcdf(netcdf_file,3)                        ! *netCDF*
        endif
C
      stop
C
      end
C
C     End of main program
C
C-----------------------------------------------------------------------
C
      subroutine advave(curv2d)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion.      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real curv2d(im,jm)
      integer i,j
C
C     u-advection and diffusion:
C
C     Advective fluxes:
C
      do j=1,jm
        do i=1,im
          advua(i,j)=0.0
        end do
      end do
C
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.1250*((d(i+1,j)+d(i,j))*ua(i+1,j)
     $                       +(d(i,j)+d(i-1,j))*ua(i,j))
     $                      *(ua(i+1,j)+ua(i,j))
        end do
      end do
C
      do j=2,jm
        do i=2,im
          fluxva(i,j)=.1250*((d(i,j)+d(i,j-1))*va(i,j)
     $                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))
     $                      *(ua(i,j)+ua(i,j-1))
        end do
      end do
C
C     Add viscous fluxes:
C
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)
     $                 -d(i,j)*2.0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))
     $                   /dx(i,j)
        end do
      end do
C
      do j=2,jm
        do i=2,im
          tps(i,j)=.250*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))
     $              *(aam2d(i,j)+aam2d(i,j-1)
     $                +aam2d(i-1,j)+aam2d(i-1,j-1))
     $              *((uab(i,j)-uab(i,j-1))
     $                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     $               +(vab(i,j)-vab(i-1,j))
     $                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.250
     $                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)
     $                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do
C
C     u-advection and diffusion:
C
      do j=1,jm
        do i=1,im
          advva(i,j)=0.0
        end do
      end do
C
C     Advective fluxes:
C
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.1250*((d(i,j)+d(i-1,j))*ua(i,j)
     $                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))
     $                      *(va(i-1,j)+va(i,j))
        end do
      end do
C
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.1250*((d(i,j+1)+d(i,j))*va(i,j+1)
     $                       +(d(i,j)+d(i,j-1))*va(i,j))
     $                      *(va(i,j+1)+va(i,j))
        end do
      end do
C
C     Add viscous fluxes:
C
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)
     $                 -d(i,j)*2.0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))
     $                   /dy(i,j)
        end do
      end do
C
      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.250
     $                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)
     $                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do
C
      if(mode.eq.2) then
C
        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.50*(cbc(i,j)+cbc(i-1,j))
     $                  *sqrt(uab(i,j)**2
     $                        +(.250*(vab(i,j)+vab(i,j+1)
     $                                 +vab(i-1,j)+vab(i-1,j+1)))**2)
     $                  *uab(i,j)
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.50*(cbc(i,j)+cbc(i,j-1))
     $                  *sqrt(vab(i,j)**2
     $                        +(.250*(uab(i,j)+uab(i+1,j)
     $                                +uab(i,j-1)+uab(i+1,j-1)))**2)
     $                  *vab(i,j)
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.250
     $                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))
     $                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))
     $                   /(dx(i,j)*dy(i,j))
          end do
        end do
C
        do j=2,jmm1
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.250
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
        end do
C
        do j=3,jmm1
          do i=2,imm1
            advva(i,j)=advva(i,j)+arv(i,j)*.250
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
        end do
C
      endif
C
      return
C
      end
C
      subroutine advct(xflux,yflux,curv)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates the horizontal portions of momentum      *
C *                advection well in advance of their use in advu and  *
C *                advv so that their vertical integrals (created in   *
C *                the main program) may be used in the external (2-D) *
C *                mode calculation.                                   *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real curv(im,jm,kb)
      real dtaam
      integer i,j,k
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            curv(i,j,k)=0.0
            advx(i,j,k)=0.0
            xflux(i,j,k)=0.0
            yflux(i,j,k)=0.0
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.250*((v(i,j+1,k)+v(i,j,k))
     $                         *(dy(i+1,j)-dy(i-1,j))
     $                         -(u(i+1,j,k)+u(i,j,k))
     $                         *(dx(i,j+1)-dx(i,j-1)))
     $                       /(dx(i,j)*dy(i,j))
          end do
        end do
      end do
C
C     Calculate x-component of velocity advection:
C
C     Calculate horizontal advective fluxes:
C
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.1250*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)
     $                           +(dt(i,j)+dt(i-1,j))*u(i,j,k))
     $                         *(u(i+1,j,k)+u(i,j,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.1250*((dt(i,j)+dt(i,j-1))*v(i,j,k)
     $                           +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))
     $                         *(u(i,j,k)+u(i,j-1,k))
          end do
        end do
      end do
C
C    Add horizontal diffusive fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.0
     $                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam=.250*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
C
            xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
            yflux(i,j,k)=.250*(dx(i,j)+dx(i-1,j)
     $                          +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do horizontal advection:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)
     $                   +yflux(i,j+1,k)-yflux(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=3,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.250
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            advy(i,j,k)=0.0
            xflux(i,j,k)=0.0
            yflux(i,j,k)=0.0
          end do
        end do
      end do
C
C     Calculate y-component of velocity advection:
C
C     Calculate horizontal advective fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.1250*((dt(i,j)+dt(i-1,j))*u(i,j,k)
     $                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))
     $                         *(v(i,j,k)+v(i-1,j,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k)=.1250*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)
     $                           +(dt(i,j)+dt(i,j-1))*v(i,j,k))
     $                         *(v(i,j+1,k)+v(i,j,k))
          end do
        end do
      end do
C
C    Add horizontal diffusive fluxes:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam=.250*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.0
     $                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)
C
            xflux(i,j,k)=.250*(dy(i,j)+dy(i-1,j)
     $                          +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do horizontal advection:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                   +yflux(i,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=3,jmm1
          do i=2,imm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.250
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advq(qb,q,qf,xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates horizontal advection and diffusion, and  *
C *                vertical advection for turbulent quantities.        *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real qb(im,jm,kb),q(im,jm,kb),qf(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
C     Do horizontal advection:
C
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.1250*(q(i,j,k)+q(i-1,j,k))
     $                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.1250*(q(i,j,k)+q(i,j-1,k))
     $                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
C
C     Do horizontal diffusion:
C
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.250*(aam(i,j,k)+aam(i-1,j,k)
     $                            +aam(i,j,k-1)+aam(i-1,j,k-1))
     $                          *(h(i,j)+h(i-1,j))
     $                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)
     $                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.250*(aam(i,j,k)+aam(i,j-1,k)
     $                            +aam(i,j,k-1)+aam(i,j-1,k-1))
     $                          *(h(i,j)+h(i,j-1))
     $                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)
     $                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.50*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.50*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
C     Do vertical advection, add flux terms, then step forward in time:
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))
     $                 *art(i,j)/(dz(k)+dz(k-1))
     $                 +xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)
     $                 *qb(i,j,k)-dti2*qf(i,j,k))
     $                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advt1(fb,f,fclim,ff,xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is centred scheme, as originally provide in    *
C *                POM (previously called advt).                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
      do j=1,jm
        do i=1,im
           f(i,j,kb)=f(i,j,kbm1)
           fb(i,j,kb)=fb(i,j,kbm1)
        end do
      end do
C
C     Do advective fluxes:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.250*((dt(i,j)+dt(i-1,j))
     $                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.250*((dt(i,j)+dt(i,j-1))
     $                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do
C
C     Add diffusive fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.50*(aam(i,j,k)+aam(i-1,j,k))
     $                         *(h(i,j)+h(i-1,j))*tprni
     $                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.50*(aam(i,j,k)+aam(i,j-1,k))
     $                         *(h(i,j)+h(i,j-1))*tprni
     $                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.50*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.50*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
C
C     Do vertical advection:
C
      do j=2,jmm1
        do i=2,imm1
          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
          zflux(i,j,kb)=0.0
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.50*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do
C
C     Add net horizontal fluxes and then step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
C
            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)
     $                 -dti2*ff(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advt2(fb,f,fclim,ff,xflux,yflux,nitera,sw)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is a first-order upstream scheme, which        *
C *                reduces implicit diffusion using the Smolarkiewicz  *
C *                iterative upstream scheme with an antidiffusive     *
C *                velocity.                                           *
C *                                                                    *
C *                It is based on the subroutines of Gianmaria Sannino *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                The calculations have been simplified by removing   *
C *                the shock switch option. It should be noted that    *
C *                this implementation does not include cross-terms    *
C *                which are in the original formulation.              *
C *                                                                    *
C *                fb,f,fclim,ff . as used in subroutine advt1         *
C *                xflux,yflux ... working arrays used to save memory  *
C *                nitera ........ number of iterations. This should   *
C *                                be in the range 1 - 4. 1 is         *
C *                                standard upstream differencing;     *
C *                                3 adds 50% CPU time to POM.         *
C *                sw ............ smoothing parameter. This should    *
C *                                preferably be 1, but 0 < sw < 1     *
C *                                gives smoother solutions with less  *
C *                                overshoot when nitera > 1.          *
C *                                                                    *
C *                Reference:                                          *
C *                                                                    *
C *                Smolarkiewicz, P.K.; A fully multidimensional       *
C *                  positive definite advection transport algorithm   *
C *                  with small implicit diffusion, Journal of         *
C *                  Computational Physics, 54, 325-362, 1984.         *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real sw
      integer nitera
      real fbmem(im,jm,kb),eta(im,jm)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      integer i,j,k,itera
C
C     Calculate horizontal mass fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            xmassflux(i,j,k)=0.0
            ymassflux(i,j,k)=0.0
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=0.250*(dy(i-1,j)+dy(i,j))
     $                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do
C
        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=0.250*(dx(i,j-1)+dx(i,j))
     $                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          fb(i,j,kb)=fb(i,j,kbm1)
          eta(i,j)=etb(i,j)
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            zwflux(i,j,k)=w(i,j,k)
            fbmem(i,j,k)=fb(i,j,k)
          end do
        end do
      end do
C
C     Start Smolarkiewicz scheme:
C
      do itera=1,nitera
C
C     Upwind advection scheme:
C
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.50
     $                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))
     $                        *fbmem(i-1,j,k)+
     $                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))
     $                        *fbmem(i,j,k))
C
              yflux(i,j,k)=0.50
     $                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j-1,k)+
     $                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j,k))
            end do
          end do
        end do
C
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,1)=0.0
            if(itera.eq.1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
            zflux(i,j,kb)=0.0
          end do
        end do
C
        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.50
     $                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k)+
     $                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do
C
C     Add net advective fluxes and step forward in time:
C
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)
     $                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
C
C     Calculate antidiffusion velocity:
C
        call smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
C
        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            do k=1,kb
              fbmem(i,j,k)=ff(i,j,k)
            end do
          end do
        end do
C
C     End of Smolarkiewicz scheme
C
      end do
C
C     Add horizontal diffusive fluxes:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.50*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.50*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni
     $                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))*0.50/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni
     $                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))*0.50/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
C
C     Add net horizontal fluxes and step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)
     $                               +yflux(i,j+1,k)-yflux(i,j,k))
     $                           /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advu
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  advu                                                *
C *                                                                    *
C * FUNCTION    :  Does horizontal and vertical advection of           *
C *                u-momentum, and includes coriolis, surface slope    *
C *                and baroclinic terms.                               *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
C     Do vertical advection:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            uf(i,j,k)=0.0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.250*(w(i,j,k)+w(i-1,j,k))
     $                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do
C
C     Combine horizontal and vertical advection with coriolis, surface
C     slope and baroclinic terms:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)
     $                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)
     $                 -aru(i,j)*.250
     $                   *(cor(i,j)*dt(i,j)
     $                      *(v(i,j+1,k)+v(i,j,k))
     $                     +cor(i-1,j)*dt(i-1,j)
     $                       *(v(i-1,j+1,k)+v(i-1,j,k)))
     $                 +grav*.1250*(dt(i,j)+dt(i-1,j))
     $                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)
     $                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.0)
     $                   *(dy(i,j)+dy(i-1,j))
     $                 +drhox(i,j,k)
          end do
        end do
      end do
C
C     Step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))
     $                 *aru(i,j)*ub(i,j,k)
     $                 -2.0*dti2*uf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))
     $                  *aru(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine advv
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Does horizontal and vertical advection of           *
C *                v-momentum, and includes coriolis, surface slope    *
C *                and baroclinic terms.                               *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
C     Do vertical advection:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            vf(i,j,k)=0.0
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.250*(w(i,j,k)+w(i,j-1,k))
     $                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
C
C     Combine horizontal and vertical advection with coriolis, surface
C     slope and baroclinic terms:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)
     $                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)
     $                 +arv(i,j)*.250
     $                   *(cor(i,j)*dt(i,j)
     $                      *(u(i+1,j,k)+u(i,j,k))
     $                     +cor(i,j-1)*dt(i,j-1)
     $                       *(u(i+1,j-1,k)+u(i,j-1,k)))
     $                 +grav*.1250*(dt(i,j)+dt(i,j-1))
     $                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)
     $                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.0)
     $                   *(dx(i,j)+dx(i,j-1))
     $                 +drhoy(i,j,k)
          end do
        end do
      end do
C
C     Step forward in time:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))
     $                 *arv(i,j)*vb(i,j,k)
     $                 -2.0*dti2*vf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     $                  *arv(i,j))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine areas_masks
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates areas and masks.                         *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j
C
C     Calculate areas of "t" and "s" cells:
C
      do j=1,jm
        do i=1,im
          art(i,j)=dx(i,j)*dy(i,j)
        end do
      end do
C
C     Calculate areas of "u" and "v" cells:
C
      do j=2,jm
        do i=2,im
          aru(i,j)=.250*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.250*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
C
      do j=1,jm
        aru(1,j)=aru(2,j)
        arv(1,j)=arv(2,j)
      end do
C
      do i=1,im
        aru(i,1)=aru(i,2)
        arv(i,1)=arv(i,2)
      end do
C
C     Initialise and set up free surface mask:
C
      do j=1,jm
        do i=1,im
          fsm(i,j)=0.0
          dum(i,j)=0.0
          dvm(i,j)=0.0
          if(h(i,j).gt.1.0) fsm(i,j)=1.0
        end do
      end do
C
C     Set up velocity masks:
C
      do j=2,jm
        do i=2,im
          dum(i,j)=fsm(i,j)*fsm(i-1,j)
          dvm(i,j)=fsm(i,j)*fsm(i,j-1)
        end do
      end do
C
      return
C
      end
C
      subroutine baropg
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j,k
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
          end do
        end do
      end do
C
C     Calculate x-component of baroclinic pressure gradient:
C
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.50*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))
     $                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                    +grav*.250*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i-1,j))
     $                      *(rho(i,j,k)-rho(i-1,j,k)
     $                        +rho(i,j,k-1)-rho(i-1,j,k-1))
     $                    +grav*.250*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i-1,j))
     $                      *(rho(i,j,k)+rho(i-1,j,k)
     $                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.250*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do
C
C     Calculate y-component of baroclinic pressure gradient:
C
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.50*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))
     $                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                    +grav*.250*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i,j-1))
     $                      *(rho(i,j,k)-rho(i,j-1,k)
     $                        +rho(i,j,k-1)-rho(i,j-1,k-1))
     $                    +grav*.250*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i,j-1))
     $                      *(rho(i,j,k)+rho(i,j-1,k)
     $                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.250*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do
C
      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine bcond(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Applies open boundary conditions.                   *
C *                                                                    *
C *                Closed boundary conditions are automatically        *
C *                enabled through specification of the masks, dum,    *
C *                dvm and fsm, in which case the open boundary        *
C *                conditions, included below, will be overwritten.    *
C *                                                                    *
C *                                The C-Grid                          *
C *                                **********                          *
C *                                                                    *
C *                The diagram is for the case where u and v are the   *
C *                primary boundary conditions together with t and     *
C *                s (co-located with el)                              * 
C *                                                                    *
C *                All interpolations are centered in space except     *
C *                those at lateral open boundary where an upstream    *
C *                Horizontal locations of e(el), t and s (etc.) are   *
C *                coincident.                                         *
C *                                                                    *
C *                People not acquainted with sigma coordinates have   *
C *                often asked what kind of boundary condition is      *
C *                applied along closed horizontal boundaries.         *
C *                Although the issue is not as important as it might  *
C *                be  for z-level grids, a direct answer is "half-    *
C *                slip" which, of course, is between free slip and    *
C *                non-slip.                                           *
C
C
C East and West end points for the C-grid in POM.
C
C                      west
C
C           v(1,j+1)=0           v(2,j+1) . . . . 
C          
C     ----<---<----<-----
C     |                 |       
C  u(1,j)   el(1,j)   u(2,j)=BC  el(2,j)   u(3,j) . . . . 
C             |                   |                 
C             -----<----<----<-----            
C
C           v(1,j)=0              v(2,j) . . . . 
C
C                                                    east
C
C                              . . . .  v(im-1,j+1)           v(im,j+1)=0
C                       
C                                           
C                 . . .  .  u(im-1,j)   el(im-1,j)  u(im,j)=BC  el(im,j)
C                                            |                   | 
C                                            ----->----->---->----
C
C                              . . . .   v(im-1,j)             v(im,j)=0
C
C  Notes:
C    1. The suffixes, f  or af, have been deleted.
C    2. All variables NOT designated as boundary condition (=BC) or set to 
C zero or obtained from an interior point are calculated points.
C    3. u(1,j) is never used but is obtained from the interior point for
C cosmetic output. Its counterpart, u(im+1,j), does not exist.
C    4. v=0 at i=1 and i=im are used as open inflow BC s unless specified
C otherwise.
C    5. The south and north extremal points are obtained from the above by 
C permuting u to v, v to u, i to j and j to i.


C **********************************************************************

      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real ga,u1,wm
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     In this example, the governing boundary conditions are a radiation
C     condition on uaf in the east and in the west, and vaf in the north
C     and south. The tangential velocities are set to zero on both
C     boundaries. These are only one set of possibilities and may not
C     represent a choice which yields the most physically realistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
C
        do j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
C
        do i=1,im
          elf(i,1)=elf(i,2)
          elf(i,jm)=elf(i,jmm1)
        end do
C
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
C
        do j=2,jmm1
C
C     East:
C
          uaf(im,j)=uabe(j)
     $               +rfe*sqrt(grav/h(imm1,j))
     $                         *(el(imm1,j)-ele(j))
          uaf(im,j)=ramp*uaf(im,j)
          vaf(im,j)=0.0
C
C     West:
C
          uaf(2,j)=uabw(j)
     $              -rfw*sqrt(grav/h(2,j))
     $                        *(el(2,j)-elw(j))
          uaf(2,j)=ramp*uaf(2,j)
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.0
C
        end do
C
        do i=2,imm1
C
C     North:
C
          vaf(i,jm)=vabn(i)
     $               +rfn*sqrt(grav/h(i,jmm1))
     $                         *(el(i,jmm1)-eln(i))
          vaf(i,jm)=ramp*vaf(i,jm)
          uaf(i,jm)=0.0
C
C     South:
C
          vaf(i,2)=vabs(i)
     $              -rfs*sqrt(grav/h(i,2))
     $                        *(el(i,2)-els(i))
          vaf(i,2)=ramp*vaf(i,2)
          vaf(i,1)=vaf(i,2)
          uaf(i,1)=0.0
C
        end do
C
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
C
C     Velocity (radiation conditions; smoothing is used in the direction
C     tangential to the boundaries):
C
        do k=1,kbm1
          do j=2,jmm1
C
C     East:
C
            ga=sqrt(h(im,j)/hmax)
            uf(im,j,k)=ga*(.250*u(imm1,j-1,k)+.50*u(imm1,j,k)
     $                     +.250*u(imm1,j+1,k))
     $                  +(1.0-ga)*(.250*u(im,j-1,k)+.50*u(im,j,k)
     $                    +.250*u(im,j+1,k))
            vf(im,j,k)=0.0
C
C     West:
C
            ga=sqrt(h(1,j)/hmax)
            uf(2,j,k)=ga*(.250*u(3,j-1,k)+.50*u(3,j,k)
     $                    +.250*u(3,j+1,k))
     $                 +(1.0-ga)*(.250*u(2,j-1,k)+.50*u(2,j,k)
     $                   +.250*u(2,j+1,k))
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.0
          end do
        end do
C
        do k=1,kbm1
          do i=2,imm1
C
C     North:
C
            ga=sqrt(h(i,jm)/hmax)
            vf(i,jm,k)=ga*(.250*v(i-1,jmm1,k)+.50*v(i,jmm1,k)
     $                     +.250*v(i+1,jmm1,k))
     $                  +(1.0-ga)*(.250*v(i-1,jm,k)+.50*v(i,jm,k)
     $                    +.250*v(i+1,jm,k))
            uf(i,jm,k)=0.0
C
C     South:
C
            ga=sqrt(h(i,1)/hmax)
            vf(i,2,k)=ga*(.250*v(i-1,3,k)+.50*v(i,3,k)
     $                    +.250*v(i+1,3,k))
     $                 +(1.0-ga)*(.250*v(i-1,2,k)+.50*v(i,2,k)
     $                   +.250*v(i+1,2,k))
            vf(i,1,k)=vf(i,2,k)
            uf(i,1,k)=0.0
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
        do k=1,kbm1
          do j=1,jm
C
C     East:
C
            u1=2.0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.0) then
              uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
            else
              uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.50*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(imm1,j))
                uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
              endif
            endif
C
C     West:
C
            u1=2.0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.0) then
              uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
            else
              uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.50*(w(2,j,k)+w(2,j,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(2,j))
                uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do i=1,im
C
C     North:
C
            u1=2.0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.0) then
              uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
            else
              uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.50*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
              endif
            endif
C
C     South:
C
            u1=2.0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.0) then
              uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
            else
              uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.50*(w(i,2,k)+w(i,2,k+1))*dti
     $              /((zz(k-1)-zz(k+1))*dt(i,2))
                uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
              endif
            endif
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
        do k=1,kb
          do j=1,jm
C
C     East:
C
            u1=2.0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.0) then
              uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
            else
              uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
            endif
C
C     West:
C
            u1=2.0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.0) then
              uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
            else
              uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
            endif
          end do
        end do
C
        do k=1,kb
          do i=1,im
C
C     North:
C
            u1=2.0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.0) then
              uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
            else
              uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
            endif
C
C     South:
C
            u1=2.0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.0) then
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
            else
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
            endif
          end do
        end do
C
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine bcondorl(idx)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  This is an optional subroutine replacing  bcond and *
C *                using Orlanski's scheme (J. Comp. Phys. 21, 251-269,*
C *                1976), specialized for the seamount problem. To     *
C *                make it work for the seamount problem, I (G.M.)     *
C *                have had to add an extra condition on an "if"       *
C *                statement in the t and s open boundary conditions,  *
C *                which involves the sign of the normal velocity.     *
C *                Thus:                                               *
C *                                                                    *
C *            if(cl.eq.0.0.and.ubw(j,k).ge.0.0) uf(1,j,k)=tbw(j,k), *
C *                                                                    *
C *                plus 3 others of the same kind.                     *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer idx
      real cl,denom
      integer i,j,k
C
      if(idx.eq.1) then
C
C-----------------------------------------------------------------------
C
C     External (2-D) boundary conditions:
C
C     In this example the governing boundary conditions are a radiation
C     condition on uaf(im,j) in the east and an inflow uaf(2,j) in the
C     west. The tangential velocities are set to zero on both
C     boundaries. These are only one set of possibilities and may not
C     represent a choice which yields the most physically realistic
C     result.
C
C     Elevation (in this application, elevation is not a primary
C     boundary condition):
C
        do  j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
C
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.2) then
C
C     External (2-D) velocity:
C
        do j=2,jmm1
C
C     West:
C
          uaf(2,j)=ramp*uabw(j)-sqrt(grav/h(2,j))*(el(2,j)-elw(j))
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.0
C
C     East:
C
          uaf(im,j)=ramp*uabe(j)
     $               +sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
          vaf(im,j)=0.0
C
        end do
C
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
C
        return
C
      else if(idx.eq.3) then
C
C-----------------------------------------------------------------------
C
C     Internal (3-D) boundary conditions:
C
C     Eastern and western radiation boundary conditions according to
C     Orlanski's explicit scheme:
C
        do k=1,kbm1
          do j=2,jmm1
C
C     West:
C
            denom=(uf(3,j,k)+ub(3,j,k)-2.0*u(4,j,k))
            if(denom.eq.0.0)denom=0.010
            cl=(ub(3,j,k)-uf(3,j,k))/denom
            if(cl.gt.1.0) cl=1.0
            if(cl.lt.0.0) cl=0.0
            uf(2,j,k)=(ub(2,j,k)*(1.0-cl)+2.0*cl*u(3,j,k))
     $                 /(1.0+cl)
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.0
C
C     East:
C
            denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.0*u(im-2,j,k))
            if(denom.eq.0.0)denom=0.010
            cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.0) cl=1.0
            if(cl.lt.0.0) cl=0.0
            uf(im,j,k)=(ub(im,j,k)*(1.0-cl)+2.0*cl*u(im-1,j,k))
     $                  /(1.0+cl)
            vf(im,j,k)=0.0
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.4) then
C
C     Temperature and salinity boundary conditions (using uf and vf,
C     respectively):
C
        do k=1,kbm1
          do j=1,jm
C
C     West:
C
            ubw(j,k)=ub(2,j,k)
            denom=(uf(2,j,k)+tb(2,j,k)-2.0*t(3,j,k))
            if(denom.eq.0.0) denom=0.010
            cl=(tb(2,j,k)-uf(2,j,k))/denom
            if(cl.gt.1.0) cl=1.0
            if(cl.lt.0.0) cl=0.0
            uf(1,j,k)=(tb(1,j,k)*(1.0-cl)+2.0*cl*t(2,j,k))/(1.0+cl)
            if(cl.eq.0.0.and.ubw(j,k).ge.0.0) uf(1,j,k)=tbw(j,k)
C
            denom=(vf(2,j,k)+sb(2,j,k)-2.0*s(3,j,k))
            if(denom.eq.0.0) denom=0.010
            cl=(sb(2,j,k)-vf(2,j,k))/denom
            if(cl.gt.1.0) cl=1.0
            if(cl.lt.0.0) cl=0.0
            vf(1,j,k)=(sb(1,j,k)*(1.0-cl)+2.0*cl*s(2,j,k))/(1.0+cl)
            if(cl.eq.0.0.and.ubw(j,k).ge.0.0) vf(1,j,k)=sbw(j,k)
C
C     East:
C
            ube(j,k)=ub(im,j,k)
            denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.0*t(im-2,j,k))
            if(denom.eq.0.0) denom=0.010
            cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.0) cl=1.0
            if(cl.lt.0.0) cl=0.0
            uf(im,j,k)=(tb(im,j,k)*(1.0-cl)+2.0*cl*t(im-1,j,k))
     $                  /(1.0+cl)
            if(cl.eq.0.0.and.ube(j,k).le.0.0) uf(im,j,k)=tbe(j,k)
C
            denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.0*s(im-2,j,k))
            if(denom.eq.0.0) denom=0.010
            cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
            if(cl.gt.1.0) cl=1.0
            if(cl.lt.0.0) cl=0.0
            vf(im,j,k)=(sb(im,j,k)*(1.0-cl)+2.0*cl*s(im-1,j,k))
     $                  /(1.0+cl)
            if(cl.eq.0.0.and.ube(j,k).le.0.0) vf(im,j,k)=sbe(j,k)
C
          end do
        end do
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.5) then
C
C     Vertical velocity boundary conditions:
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      else if(idx.eq.6) then
C
C     q2 and q2l boundary conditions:
C
        do k=1,kb
C
          do j=1,jm
            uf(im,j,k)=1.e-10
            vf(im,j,k)=1.e-10
            uf(1,j,k)=1.e-10
            vf(1,j,k)=1.e-10
          end do
C
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
C
        return
C
      endif
C
      end
C
      subroutine box
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up conservation box problem.                   *
C *                                                                    *
C *                This basin uses the same grid as the seamount       *
C *                problem, but it has a flat bottom, is surrounded by *
C *                walls and is initialised with uniform salinity and  *
C *                temperature. It is forced by a surface input of     *
C *                water of the same temperature and salinity as the   *
C *                water in the basin. Therefore, the temperature and  *
C *                salinity in the basin should not change, and the    *
C *                free surface should fall at a rate vflux. It is also*
C *                forced by a steady atmospheric pressure field which *
C *                depresses the southwestern half of the model by 1 m *
C *                and elevates the northeastern half of the model by  *
C *                1 m.                                                *
C *                                                                    *
C *                Since this problem defines its own fixed e_atmos,   *
C *                tatm, satm and e_atmos, comment out corresponding   *
C *                declarations after the do 9000 statement in main    *
C *                program.                                            *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real depth,delx,tatm,satm
      integer i,j,k
C
C     Water depth:
C
      depth=4500.0
C
C     Grid size:
C
      delx=8000.0
C
C     Set up grid dimensions, areas of free surface cells, and
C     Coriolis parameter:
C
      do j=1,jm
        do i=1,im
C
C     For constant grid size:
C
C         dx(i,j)=delx
C         dy(i,j)=delx
C
C     For variable grid size:
C
          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.0
          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.0
C
          cor(i,j)=1.e-4
C
        end do
      end do
C
C     Calculate horizontal coordinates of grid points and rotation
C     angle.
C
C     NOTE that this is introduced solely for the benefit of any post-
C     processing software, and in order to conform with the requirements
C     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
C
C     There are four horizontal coordinate systems, denoted by the
C     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
C     "e" is an elevation point and "c" is a cell corner), as shown
C     below. In addition, "east_*" is an easting and "north_*" is a
C     northing. Hence the coordinates of the "u" points are given by
C     (east_u,north_u).
C
C     Also, if the centre point of the cell shown below is at
C     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
C     the coordinates of the western of the two "u" points,
C     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
C     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
C     coordinates of the southwestern corner point of the cell. The
C     southwest corner of the entire grid is at
C     (east_c(1,1),north_c(1,1)).
C
C                      |              |
C                    --c------v-------c--
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                      u      e       u
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                    --c------v-------c--
C                      |              |
C
C
C     NOTE that the following calculation of east_c and north_c only
C     works properly for a rectangular grid with east and north aligned
C     with i and j, respectively:
C
      do j=1,jm
        east_c(1,j)=0.0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
C
      do i=1,im
        north_c(i,1)=0.0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
C
C     The following works properly for any grid:
C
C     Elevation points:
C
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)
     $                  +east_c(i,j+1)+east_c(i+1,j+1))/4.0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)
     $                   +north_c(i,j+1)+north_c(i+1,j+1))/4.0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im-1
        east_e(i,jm)
     $    =((east_c(i,jm)+east_c(i+1,jm))*3.0
     $       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.0
        north_e(i,jm)
     $    =((north_c(i,jm)+north_c(i+1,jm))*3.0
     $       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.0
      end do
C
      do j=1,jm-1
        east_e(im,j)
     $    =((east_c(im,j)+east_c(im,j+1))*3.0
     $       -east_c(im-1,j)-east_c(im-1,j+1))/4.0
        north_e(im,j)
     $    =((north_c(im,j)+north_c(im,j+1))*3.0
     $       -north_c(im-1,j)-north_c(im-1,j+1))/4.0
      end do
C
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)
     $               -(east_e(im-2,jm)+east_e(im,jm-2))/2.0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)
     $               -(north_e(im-2,jm)+north_e(im,jm-2))/2.0
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.0-east_c(i,jm-1))/2.0
        north_u(i,jm)=(north_c(i,jm)*3.0-north_c(i,jm-1))/2.0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.0-east_c(im-1,j))/2.0
        north_v(im,j)=(north_c(im,j)*3.0-north_c(im-1,j))/2.0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre:
C
C     (NOTE that the following calculation of rot only works properly
C     for this particular rectangular grid)
C
      do j=1,jm
        do i=1,im
          rot(i,j)=0.0
        end do
      end do
C
C     Define depth:
C
      do i=1,im
        do j=1,jm
          h(i,j)=depth
        end do
      end do
C
C     Close the north and south boundaries:
C
      do i=1,im
        h(i,1)=1.0
        h(i,jm)=1.0
      end do
C
C     Close the east and west boundaries:
C
      do j=1,jm
        h(1,j)=1.0
        h(im,j)=1.0
      end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
      if(slmax.lt.1.0) call slpmax
C
C     Set tbias and sbias here for test (tbias and sbias would
C     normally only be set in the main program):
C
      tbias=10.0
      sbias=20.0
      write(6,1) tbias,sbias
    1 format(/' tbias and sbias changed in subroutine box to:'/
     $         2f10.3//)
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=20.0-tbias
            sb(i,j,k)=35.0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
          end do
        end do
      end do
C
C     Initialise uab and vab as necessary
C     (NOTE that these have already been initialised to zero in the
C     main program):
C
      do j=1,jm
        do i=1,im
C     No conditions necessary for this problem
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
      do j=1,jm
        do i=1,im
          if(i+j-57.le.0) then
            e_atmos(i,j)=1.0
          else
            e_atmos(i,j)=-1.0
          endif
C
C     Ensure atmospheric pressure cannot make water depth go negative:
C
          e_atmos(i,j)=min(e_atmos(i,j),h(i,j))
C
          vfluxf(i,j)=-0.00010
C
C     See main program, just after "Begin numerical integration", for
C     an explanation of these terms:
C 
          tatm=20.0
          satm=35.0
C
        end do
      end do
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in this problem, all lateral boundaries are closed through
C     the specification of the masks fsm, dum and dvm):
C
      rfe=1.0
      rfw=1.0
      rfn=1.0
      rfs=1.0
C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
C
      end
C
      subroutine dens(si,ti,rhoo)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates (density-1000.)/rhoref.                  *
C *                                                                    *
C *                (see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech.,  *
C *                609-611.)                                           *
C *                                                                    *
C *                ti is potential temperature                         *
C *                                                                    *
C *                If using 32 bit precision, it is recommended that   *
C *                cr,p,rhor,sr,tr,tr2,tr3 and tr4 be made double      *
C *                precision, and the "e"s in the constants be changed *
C *                to "d"s.                                            *
C *                                                                    *
C * NOTE: if pressure is not used in dens, buoyancy term (boygr)       *
C *       in profq must be changed (see note in profq)                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real si(im,jm,kb),ti(im,jm,kb),rhoo(im,jm,kb)
      real cr,p,rhor,sr,tr,tr2,tr3,tr4
      integer i,j,k
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
C
            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr
C
C     Approximate pressure in units of bars:
C
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5
C
            rhor=-0.1574060+6.793952e-2*tr
     $            -9.095290e-3*tr2+1.001685e-4*tr3
     $            -1.120083e-6*tr4+6.536332e-9*tr4*tr
C
            rhor=rhor+(0.8244930-4.0899e-3*tr
     $            +7.6438e-5*tr2-8.2467e-7*tr3
     $            +5.3875e-9*tr4)*sr
     $            +(-5.72466e-3+1.0227e-4*tr
     $            -1.6546e-6*tr2)*abs(sr)**1.5
     $            +4.8314e-4*sr*sr
C
            cr=1449.10+.08210*p+4.550*tr-.0450*tr2
     $          +1.340*(sr-35.0)
            rhor=rhor+1.e5*p/(cr*cr)*(1.0-2.0*p/(cr*cr))
C
            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)
C
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine depth
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Establishes the vertical sigma grid with log        *
C *                distributions at the top and bottom and a linear    *
C *                distribution in between. The number of layers of    *
C *                reduced thickness are kl1-2 at the surface and      *
C *                kb-kl2-1 at the bottom. kl1 and kl2 are defined in  *
C *                the main program. For no log portions, set kl1=2    *
C *                and kl2=kb-1.                                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real delz
      integer kdz(12)
      integer k
C
      data kdz/1,1,2,4,8,16,32,64,128,256,512,1024/
C
      z(1)=0.0
C
      do k=2,kl1
        z(k)=z(k-1)+kdz(k-1)
      end do
C
      delz=z(kl1)-z(kl1-1)
C
      do k=kl1+1,kl2
        z(k)=z(k-1)+delz
      end do
C
      do k=kl2+1,kb
        dz(k)=float(kdz(kb-k+1))*delz/float(kdz(kb-kl2))
        z(k)=z(k-1)+dz(k)
      end do
C
      do k=1,kb
        z(k)=-z(k)/z(kb)
      end do
C
      do k=1,kb-1
        zz(k)=0.50*(z(k)+z(k+1))
      end do
C
      zz(kb)=2.0*zz(kb-1)-zz(kb-2)
C
      do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
C
      dz(kb)=0.0
      dzz(kb)=0.0
C
      write(6,1)
    1 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
C
      do k=1,kb
        write(6,2) k,z(k),zz(k),dz(k),dzz(k)
    2   format((' ',i5,4f10.5))
      end do
C
      write(6,3)
    3 format(//)
C
      return
C
      end
C
      subroutine findpsi
C **********************************************************************
C *                                                                    *
C * ROUTINE NAME:  findpsi                                             *
C *                                                                    *
C * FUNCTION    :  Calculates the stream function, first assuming      *
C *                zero on the southern boundary and then, using the   *
C *                values on the western boundary, the stream function *
C *                is calculated again. If the elevation field is near *
C *                steady state, the two calculations should agree;    *
C *                otherwise not.                                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      integer i,j
C
      do j=1,jm
        do i=1,im
          psi(i,j)=0.0
        end do
      end do
C
C     Sweep northward:
C
      do j=2,jmm1
        do i=2,im
          psi(i,j+1)=psi(i,j)
     $                +.250*uab(i,j)*(d(i,j)+d(i-1,j))
     $                  *(dy(i,j)+dy(i-1,j))
        end do
      end do
C
      call prxy('Streamfunction, psi from u              ',
     $          time,psi,im,iskp,jm,jskp,0.0)
C
C    Sweep eastward:
C
      do j=2,jm
        do i=2,imm1
          psi(i+1,j)=psi(i,j)
     $                -.250*vab(i,j)*(d(i,j)+d(i,j-1))
     $                  *(dx(i,j)+dx(i,j-1))
        end do
      end do
C
      call prxy('Streamfunction, psi from v              ',
     $          time,psi,im,iskp,jm,jskp,0.0)
C
      return
C
      end
C
      subroutine file2ic
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up my own problem.                             *
C *                                                                    *
C * This example read IC from IC.dat file, generated by GRID.f in      *
C * GRID-DATA directory. Only minimal number of fields are read,       *
C * while others are calculated here.                                  *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real rad,re,dlat,dlon,cff
      integer i,j,k,m
      character*5 field
      rad=0.01745329
      re=6371.E3
C
      write(6,'(/,'' Read grid and initial conditions '',/)')
C
C--- 1D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') z
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') zz
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') dz
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') dzz
C--- 2D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') east_e
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') north_e
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') h
C--- 3D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') t
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') s
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') rmean
C--- Constant wind stress read here
C (for time dep. read in loop 9000 & interpolate in time)
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') wusurf
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') wvsurf
C
C --- print vertical grid distribution
C
      write(6,2)
    2 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
      write(6,'(''  '',/)')
      do k=1,kb
        write(6,3) k,z(k),zz(k),dz(k),dzz(k)
    3   format((' ',i5,4f10.3))
      end do
      write(6,'(''  '',//)')
C
C --- calc. surface & lateral BC from climatology
C
        do j=1,jm
          do i=1,im
             tsurf(i,j)=t(i,j,1)
             ssurf(i,j)=s(i,j,1)
            do k=1,kb
              tclim(i,j,k)=t(i,j,k)
              sclim(i,j,k)=s(i,j,k)
            end do
          end do
        end do
C
C                    --- EAST & WEST BCs ---
        do j=1,jm
              ele(j)=0.
              elw(j)=0.
C --- other vel. BCs (fixed in time) can be specified here
              uabe(j)=0.
              uabw(j)=0.
            do k=1,kb
              ubw(j,k)=0.
              ube(j,k)=0.
              tbw(j,k)=tclim(1,j,k)
              sbw(j,k)=sclim(1,j,k)
              tbe(j,k)=tclim(im,j,k)
              sbe(j,k)=sclim(im,j,k)
            end do
        end do
C                    --- NORTH & SOUTH BCs ---
        do i=1,im
              els(i)=0.
              eln(i)=0.
              vabs(i)=0.
              vabn(i)=0.
            do k=1,kb
              vbs(i,k)=0.
              vbn(i,k)=0.
              tbs(i,k)=tclim(i,1,k)
              sbs(i,k)=sclim(i,1,k)
              tbn(i,k)=tclim(i,jm,k)
              sbn(i,k)=sclim(i,jm,k)
            end do
        end do
C
C     Set initial conditions:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            tb(i,j,k)=t(i,j,k)
            sb(i,j,k)=s(i,j,k)
            ub(i,j,k)=0.
            vb(i,j,k)=0.
          end do
        end do
      end do
C
      call dens(sb,tb,rho)
C
C --- calc. Curiolis Parameter
C
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
            aam2d(i,j)=aam(i,j,1)
            elb(i,j)=0.
            etb(i,j)=0.
            dt(i,j)=h(i,j)
          end do
        end do
C
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
C
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
C
C --- the following grids are needed only for netcdf plotting
C
C     Corner of cell points:
C
      do j=2,jm
        do i=2,im
          east_c(i,j)=(east_e(i,j)+east_e(i-1,j)
     $                  +east_e(i,j-1)+east_e(i-1,j-1))/4.0
          north_c(i,j)=(north_e(i,j)+north_e(i-1,j)
     $                   +north_e(i,j-1)+north_e(i-1,j-1))/4.0
        end do
      end do
C
C
C     Extrapolate ends (approx.):
C
      do i=2,im
        east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
        north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
      end do
        east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
C
      do j=2,jm
        east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
        north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
      end do
        north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.0-east_c(i,jm-1))/2.0
        north_u(i,jm)=(north_c(i,jm)*3.0-north_c(i,jm-1))/2.0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.0-east_c(im-1,j))/2.0
        north_v(im,j)=(north_c(im,j)*3.0-north_c(im-1,j))/2.0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre: (only needed for CDF plotting)
C
      do j=1,jm
        do i=1,im-1
          rot(i,j)=0.
          dlat=north_e(i+1,j)-north_e(i,j)
          dlon= east_e(i+1,j)- east_e(i,j)
           if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
        end do
       rot(im,j)=rot(im-1,j)
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     set all=0 for closed BCs.
C     Values=0 for vel BC only, =1 is combination of vel+elev.
      rfe=0.0
      rfw=0.0
      rfn=0.0
      rfs=0.0
C
      return
      end
C
      subroutine printall
C **********************************************************************
C *                                                                    *
C *                         POM2K SOURCE CODE                          *
C *                                                                    *
C * ROUTINE NAME:  printall                                            *
C *                                                                    *
C * FUNCTION    :  Prints a set of outputs to device 6                 *
C *                                                                    *
C *                Edit as approriate.                                 *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer io(100),jo(100),ko(100)
C
      include 'pom2k.c'
C
C     2-D horizontal fields:
C
          call prxy('Depth-averaged u, uab                   ',
     $              time,uab,im,iskp,jm,jskp,0.0)
C
          call prxy('Depth-averaged v, vab                   ',
     $              time,vab,im,iskp,jm,jskp,0.0)
C
          call prxy('Surface elevation, elb                  ',
     $              time,elb,im,iskp,jm,jskp,0.0)
C
c         call prxy(' egf ',time,egf,im,iskp,jm,jskp,0.0)
c         call prxy(' utf ',time,utf,im,iskp,jm,jskp,0.0)
c         call prxy(' vtf ',time,vtf,im,iskp,jm,jskp,0.0)
c
C     Calculate and print streamfunction:
C
          call findpsi
C
          if(mode.ne.2) then
C
C     2-D horizontal sections of 3-D fields:
C
C     Set levels for output:
C
            ko(1)=1
            ko(2)=kb/2
            ko(3)=kb-1
C
            call prxyz('x-velocity, u                           ',
     $                 time,u    ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
            call prxyz('y-velocity, v                           ',
     $                 time,v    ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
            ko(1)=2
            call prxyz('z-velocity, w                           ',
     $                 time,w    ,im,iskp,jm,jskp,kb,ko,3,0.0 )
            ko(1)=1
C
            call prxyz('Potential temperature, t                ',
     $                 time,t    ,im,iskp,jm,jskp,kb,ko,3,1.e-2)
C
            call prxyz('Salinity, s                              ',
     $                 time,s    ,im,iskp,jm,jskp,kb,ko,3,1.e-2)
C
            call prxyz('(density-1000)/rhoref, rho              ',
     $                 time,rho  ,im,iskp,jm,jskp,kb,ko,3,1.e-5)
C
c           call prxyz('Turbulent kinetic energy x 2, q2        ',
c    $                 time,q2   ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
c           call prxyz('Turbulent length scale, l               ',
c    $                 time,l    ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
            call prxyz('Horizontal kinematic viscosity, aam     ',
     $                 time,aam  ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
            call prxyz('Vertical kinematic viscosity, km        ',
     $                 time,km   ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
c           call prxyz('Vertical kinematic diffusivity, kh      ',
c    $                 time,kh   ,im,iskp,jm,jskp,kb,ko,3,0.0 )
C
C     Vertical sections of 3-D fields, normal to j-axis:
C
C     Set sections for output:
C
            jo(1)=1
            jo(2)=jm/2
            jo(3)=jm-1
C
            call prxz('x-velocity, u                           ',
     $                time,u    ,im,iskp,jm,kb,jo,3,0.0 ,dt,zz)
C
            call prxz('y-velocity, v                           ',
     $                time,v    ,im,iskp,jm,kb,jo,3,0.0 ,dt,zz)
C
            call prxz('z-velocity, w                           ',
     $                time,w    ,im,iskp,jm,kb,jo,3,0.0 ,dt,z )
C
            call prxz('Potential temperature, t                ',
     $                time,t    ,im,iskp,jm,kb,jo,3,1.e-2,dt,zz)
C
            call prxz('Salinity, s                             ',
     $                time,s    ,im,iskp,jm,kb,jo,3,1.e-2,dt,zz)
C
            call prxz('(density-1000)/rhoref, rho              ',
     $                time,rho  ,im,iskp,jm,kb,jo,3,1.e-5,dt,zz)
C
c           call prxz('Turbulent kinetic energy x 2, q2        ',
c    $                time,q2   ,im,iskp,jm,kb,jo,3,0.0 ,dt,z )
C
c           call prxz('Turbulent length scale, l               ',
c    $                time,l    ,im,iskp,jm,kb,jo,3,0.0 ,dt,z )
C
c           call prxz('Horizontal kinematic viscosity, aam     ',
c    $                time,aam  ,im,iskp,jm,kb,jo,3,0.0 ,dt,zz)
C
c           call prxz('Vertical kinematic viscosity, km        ',
c    $                time,km   ,im,iskp,jm,kb,jo,3,0.0 ,dt,z )
C
c           call prxz('Vertical kinematic diffusivity, kh      ',
c    $                time,kh   ,im,iskp,jm,kb,jo,3,0.0 ,dt,z )
C
C     Vertical sections of 3-D fields, normal to i-axis:
C
C     Set sections for output:
C
            io(1)=1
            io(2)=im/2
            io(3)=im-1
C
            call pryz('x-velocity, u                           ',
     $                time,u    ,im,jm,jskp,kb,io,3,0.0 ,dt,zz)
C
            call pryz('y-velocity, v                           ',
     $                time,v    ,im,jm,jskp,kb,io,3,0.0 ,dt,zz)
C
            call pryz('z-velocity, w                           ',
     $                time,w    ,im,jm,jskp,kb,io,3,0.0 ,dt,zz)
C
            call pryz('Potential temperature, t                ',
     $                time,t    ,im,jm,jskp,kb,io,3,1.e-2,dt,zz)
C
c           call pryz('Salinity x rho / rhoref, s              ',
c    $                time,s    ,im,jm,jskp,kb,io,3,1.e-2,dt,zz)
C
c           call pryz('(density-1000)/rhoref, rho              ',
c    $                time,rho  ,im,jm,jskp,kb,io,3,1.e-5,dt,zz)
C
c           call pryz('Turbulent kinetic energy x 2, q2        ',
c    $                time,q2   ,im,jm,jskp,kb,io,3,0.0 ,dt,z )
C
c           call pryz('Turbulent length scale, l               ',
c    $                time,l    ,im,jm,jskp,kb,io,3,0.0 ,dt,z )
C
c           call pryz('Horizontal kinematic viscosity, aam     ',
c    $                time,aam  ,im,jm,jskp,kb,io,3,0.0 ,dt,zz)
C
c           call pryz('Vertical kinematic viscosity, km        ',
c    $                time,km   ,im,jm,jskp,kb,io,3,0.0 ,dt,z )
C
c           call pryz('Vertical kinematic diffusivity, kh      ',
c    $                time,kh   ,im,jm,jskp,kb,io,3,0.0 ,dt,z )
C
          endif
C
      return
C
      end
C
      subroutine profq(sm,sh,dh,cc)
C **********************************************************************
C *                                        Updated: Sep. 24, 2003      *
C * FUNCTION    :  Solves for q2 (twice the turbulent kinetic energy), *
C *                q2l (q2 x turbulent length scale), km (vertical     *
C *                kinematic viscosity) and kh (vertical kinematic     *
C *                diffusivity), using a simplified version of the     *
C *                level 2 1/2 model of Mellor and Yamada (1982).      *
C * In this version, the Craig-Banner sub-model whereby breaking wave  * 
C * tke is injected into the surface is included. However, we use an   *
C * analytical solution to the near surface tke equation to solve for  *
C * q2 at the surface giving the same result as C-B diffusion. The new *
C * scheme is simpler and more robust than the latter scheme.          *     
C *                                                                    *
C * References                                                         *
C *   Craig, P. D. and M. L. Banner, Modeling wave-enhanced turbulence *
C *     in the ocean surface layer. J. Phys. Oceanogr., 24, 2546-2559, *
C *     1994.                                                          *
C *   Ezer, T., On the seasonal mixed-layer simulated by a basin-scale *
C *     ocean model and the Mellor-Yamada turbulence scheme,           *
C *     J. Geophys. Res., 105(C7), 16,843-16,855, 2000.                *
C *   Mellor, G.L. and T. Yamada, Development of a turbulence          *
C *     closure model for geophysical fluid fluid problems,            *
C *     Rev. Geophys. Space Phys., 20, 851-875, 1982.                  *
C *   Mellor, G. L., One-dimensional, ocean surface layer modeling,    *
C *     a problem and a solution. J. Phys. Oceanogr., 31(3), 790-809,  *
C *     2001.                                                          *
C *   Mellor, G.L. and A. Blumberg, Wave breaking and ocean surface    *
C *     thermal response, J. Phys. Oceanogr., 2003.                    *
C *   Stacey, M. W., Simulations of the wind-forced near-surface       *
C *     circulation in Knight Inlet: a parameterization of the         *
C *     roughness length. J. Phys. Oceanogr., 29, 1363-1367, 1999.     *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real sm(im,jm,kb),sh(im,jm,kb),cc(im,jm,kb)
      real gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
      real prod(im,jm,kb),kn(im,jm,kb)
      real a1,a2,b1,b2,c1
      real coef1,coef2,coef3,coef4,coef5
      real const1,e1,e2,ghc
      real p,sef,sp,tp
      real l0(im,jm)
      real cbcnst,surfl,shiw
      real utau2, df0,df1,df2 
C
      integer i,j,k,ki
C
      equivalence (prod,kn)
C
      data a1,b1,a2,b2,c1/0.920,16.60,0.740,10.10,0.080/
      data e1/1.80/,e2/1.330/
      data sef/1.0/
      data cbcnst/100./surfl/2.e5/shiw/0.0/
C
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.0*umol)*.50
     $                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.0*umol)*.50
     $                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
C
C     Surface and bottom boundary conditions:
C
      const1=(16.60**(2.0/3.0))*sef
C
C initialize fields that are not calculated on all boundaries
C but are later used there
      do i=1,im
        ee(i,jm,1)=0.
        gg(i,jm,1)=0.
        l0(i,jm)=0.
      end do
      do j=1,jm
        ee(im,j,1)=0.
        gg(im,j,1)=0.
        l0(im,j)=0.
      end do
      do i=1,im
      do j=1,jm
       do k=2,kbm1
        prod(i,j,k)=0.
       end do
      end do
      end do
C
      do j=1,jmm1
        do i=1,imm1
          utau2=sqrt((.50*(wusurf(i,j)+wusurf(i+1,j)))**2
     $                  +(.50*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
C Wave breaking energy- a variant of Craig & Banner (1994)
C see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.0
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2 
C Surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2/grav
C
          uf(i,j,kb)=sqrt((.50*(wubot(i,j)+wubot(i+1,j)))**2
     $                   +(.50*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
C
C    Calculate speed of sound squared:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
C
C     Calculate pressure in units of decibars:
C
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-4
            cc(i,j,k)=1449.10+.008210*p+4.550*tp -.0450*tp**2
     $                 +1.340*(sp-35.00)
            cc(i,j,k)=cc(i,j,k)
     $                 /sqrt((1.0-.016420*p/cc(i,j,k))
     $                   *(1.0-0.400*p/cc(i,j,k)**2))
          end do
        end do
      end do
C
C     Calculate buoyancy gradient:
C
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $                    /(dzz(k-1)* h(i,j))
C *** NOTE: comment out next line if dens does not include pressure
     $      +(grav**2)*2.0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
            gh(i,j,k)=min(gh(i,j,k),.0280)
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.0
          gh(i,j,1)=0.0
          gh(i,j,kb)=0.0
        end do
      end do
C
C    Calculate production of turbulent kinetic energy:
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.250*sef
     $                   *((u(i,j,k)-u(i,j,k-1)
     $                      +u(i+1,j,k)-u(i+1,j,k-1))**2
     $                     +(v(i,j,k)-v(i,j,k-1)
     $                      +v(i,j+1,k)-v(i,j+1,k-1))**2)
     $                   /(dzz(k-1)*dh(i,j))**2
C   Add shear due to internal wave field
     $             -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do
C
C  NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000),
C  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
      ghc=-6.00
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.0
C It is unclear yet if diss. corr. is needed when surf. waves are included.
c           if(gh(i,j,k).lt.0.0)
c    $        stf(i,j,k)=1.00-0.90*(gh(i,j,k)/ghc)**1.50
c           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.10
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                   /(b1*l(i,j,k)+small)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.0/(a(i,j,k)+c(i,j,k)*(1.0-ee(i,j,k-1))
     $                      -(2.0*dti2*dtef(i,j,k)+1.0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $                 -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
C
      do j=1,jm
        do i=1,im
          ee(i,j,2)=0.0
          gg(i,j,2)=0.0
          vf(i,j,kb)=0.0
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)
     $                   *(1.0+e2*((1.0/abs(z(k)-z(1))
     $                               +1.0/abs(z(k)-z(kb)))
     $                                *l(i,j,k)/(dh(i,j)*kappa))**2)
            gg(i,j,k)=1.0/(a(i,j,k)+c(i,j,k)*(1.0-ee(i,j,k-1))
     $                      -(dti2*dtef(i,j,k)+1.0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
              uf(i,j,k)=small
              vf(i,j,k)=0.1*dt(i,j)*small
            endif
          end do
        end do
      end do
C
C-----------------------------------------------------------------------
C
C     The following section solves for km and kh:
C
      coef4=18.0*a1*a1+9.0*a1*a2
      coef5=9.0*a1*a2
C
C     Note that sm and sh limit to infinity when gh approaches 0.0288:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.0-6.0*a1/b1*stf(i,j,k))
            coef2=3.0*a2*b2/stf(i,j,k)+18.0*a1*a2
            coef3=a1*(1.0-3.0*c1-6.0*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.0-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.0-coef5*gh(i,j,k))
          end do
        end do
      end do
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            kn(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kq(i,j,k)=(kn(i,j,k)*.410*sh(i,j,k)+kq(i,j,k))*.50
            km(i,j,k)=(kn(i,j,k)*sm(i,j,k)+km(i,j,k))*.50
            kh(i,j,k)=(kn(i,j,k)*sh(i,j,k)+kh(i,j,k))*.50
          end do
        end do
      end do
C cosmetics: make boundr. values as interior
C (even if not used, printout otherwise may show strange values)
      do k=1,kb
        do i=1,im
           km(i,jm,k)=km(i,jmm1,k)*fsm(i,jm)
           kh(i,jm,k)=kh(i,jmm1,k)*fsm(i,jm)
           km(i,1,k)=km(i,2,k)*fsm(i,1)
           kh(i,1,k)=kh(i,2,k)*fsm(i,1)
        end do
        do j=1,jm
           km(im,j,k)=km(imm1,j,k)*fsm(im,j)
           kh(im,j,k)=kh(imm1,j,k)*fsm(im,j)
           km(1,j,k)=km(2,j,k)*fsm(1,j)
           kh(1,j,k)=kh(2,j,k)*fsm(1,j)
        end do
      end do
C
      return
C
      end
C
c ---------------------------------------------------------------------
C
      subroutine proft(f,wfsurf,fsurf,nbc,dh)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Solves for vertical diffusion of temperature and    *
C *                salinity using method described by Richmeyer and    *
C *                Morton.                                             *
C *                                                                    *
C *                Irradiance parameters are from Paulson and Simpson. *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                Paulson, C. A., and J. Simpson, 1977: Irradiance    *
C *                  measurements in the upper ocean, J. Phys.         *
C *                  Oceanogr., 7, 952-956.                            *
C *                                                                    *
C *                NOTES:                                              *
C *                                                                    *
C *                (1) wfsurf and swrad are negative values when water *
C *                    column is warming or salt is being added.       *
C *                                                                    *
C *                (2) nbc may only be 1 and 3 for salinity.           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real f(im,jm,kb),wfsurf(im,jm)
      real fsurf(im,jm),dh(im,jm)
      integer nbc
      real rad(im,jm,kb),r(5),ad1(5),ad2(5)
      integer i,j,k,ki
C
C-----------------------------------------------------------------------
C
C     Irradiance parameters after Paulson and Simpson:
C
C       ntp               1      2       3       4       5
C   Jerlov type           i      ia      ib      ii     iii
C
      data r   /       .580,  .620,  .670,  .770,  .780 /
      data ad1 /       .350,  .600,  1.00,  1.50,  1.40 /
      data ad2 /       23.0,  20.0,  17.0,  14.0,  7.90 /
C
C-----------------------------------------------------------------------
C
C     Surface boundary condition:
C
C       nbc   prescribed    prescribed   short wave
C             temperature      flux      penetration
C             or salinity               (temperature
C                                           only)
C
C        1        no           yes           no
C        2        no           yes           yes
C        3        yes          no            no
C        4        yes          no            yes
C
C     NOTE that only 1 and 3 are allowed for salinity.
C
C-----------------------------------------------------------------------
C
C     The following section solves the equation:
C
C       dti2*(kh*f')'-f=-fb
C
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
C     Calculate penetrative radiation. At the bottom any unattenuated
C     radiation is deposited in the bottom layer:
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            rad(i,j,k)=0.0
          end do
        end do
      end do
C
      if(nbc.eq.2.or.nbc.eq.4) then
C
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=swrad(i,j)
     $                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp))
     $                      +(1.0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
            end do
          end do
        end do
C
      endif
C
      if(nbc.eq.1) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.0)
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.0)
          end do
        end do
C
      else if(nbc.eq.2) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.0)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))
     $                 /(dz(1)*dh(i,j))
     $                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.0)
          end do
        end do
C
      else if(nbc.eq.3.or.nbc.eq.4) then
C
        do j=1,jm
          do i=1,im
            ee(i,j,1)=0.0
            gg(i,j,1)=fsurf(i,j)
          end do
        end do
C
      endif
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.0/(a(i,j,k)+c(i,j,k)*(1.0-ee(i,j,k-1))-1.0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)
     $                 +dti2*(rad(i,j,k)-rad(i,j,k+1))
     $                   /(dh(i,j)*dz(k)))
     $                 *gg(i,j,k)
          end do
        end do
      end do
C
C     Bottom adiabatic boundary condition:
C
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)
     $                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))
     $                   /(dh(i,j)*dz(kbm1)))
     $                 /(c(i,j,kbm1)*(1.0-ee(i,j,kbm2))-1.0)
        end do
      end do
C
      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
          f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine profu
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Solves for vertical diffusion of x-momentum using   *
C *                method described by Richmeyer and Morton.           *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                NOTE that wusurf has the opposite sign to the wind  *
C *                speed.                                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
      real dh(im,jm)
      integer i,j,k,ki
C
C     The following section solves the equation:
C
C       dti2*(km*u')'-u=-ub
C
      do j=1,jm
        do i=1,im
          dh(i,j)=1.0
        end do
      end do
C
      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.50
        end do
      end do
C
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.50
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.0)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))
     $               -uf(i,j,1))
     $               /(a(i,j,1)-1.0)
        end do
      end do
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.0/(a(i,j,k)+c(i,j,k)*(1.0-ee(i,j,k-1))-1.0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.50*(cbc(i,j)+cbc(i-1,j))
     $              *sqrt(ub(i,j,kbm1)**2
     $                +(.250*(vb(i,j,kbm1)+vb(i,j+1,kbm1)
     $                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.0
     $                    -(ee(i,j,kbm2)-1.0)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do
C

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
C
      return
C
      end
C
      subroutine profv
C **********************************************************************
C                                                                      *
C * FUNCTION    :  Solves for vertical diffusion of y-momentum using   *
C *                method described by Richmeyer and Morton.           *
C *                                                                    *
C *                See:                                                *
C *                                                                    *
C *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
C *                  Methods for Initial-Value Problems, 2nd edition,  *
C *                  Interscience, New York, 198-201.                  *
C *                                                                    *
C *                NOTE that wvsurf has the opposite sign to the wind  *
C *                speed.                                              *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
      real dh(im,jm)
      integer i,j,k,ki
C
C     The following section solves the equation:
C
C       dti2*(km*u')'-u=-ub
C
      do j=1,jm
        do i=1,im
          dh(i,j)=1.0
        end do
      end do
C
      do j=2,jm
        do i=2,im
          dh(i,j)=.50*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do
C
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.50
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
C
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.0)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))
     $               /(a(i,j,1)-1.0)
        end do
      end do
C
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.0/(a(i,j,k)+c(i,j,k)*(1.0-ee(i,j,k-1))-1.0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.50*(cbc(i,j)+cbc(i,j-1))
     $              *sqrt((.250*(ub(i,j,kbm1)+ub(i+1,j,kbm1)
     $                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2
     $                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.0
     $                    -(ee(i,j,kbm2)-1.0)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do
C
      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do
C
      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
C
      return
C
      end
C
      subroutine prxy(label,time,a,im,iskp,jm,jskp,scala)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes a horizontal 2-D field.                      *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jskp ........ skipping interval for j               *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm
      real a(im,jm)
      real time,scala
      integer iskp,jskp
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,jwr,cols
C
      if(scala.ge.0.0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.0) scale = 1.0
      if (scala.eq.0.0) then
        amx=1.e-12
        do j=1,jm,jskp
          do i=1,im,iskp
            amx=max(abs(a(i,j)),amx)
          end do
        end do
          scale=10.0**(int(log10(amx)+100.0)-103)
        endif
      if(scala.gt.0.0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do ib=1,im,cols*iskp
C
        ie=ib+(cols-1)*iskp
        if(ie.gt.im) ie=im
C
        if(scala.ge.0.0) then
          write(6,3) (i,i=ib,ie,iskp)
    3     format(/,2x,24i5,/)
        else
          write(6,4) (i,i=ib,ie,iskp)
    4     format(/,12i10,/)
        endif
C
        do j=1,jm,jskp
          jwr=jm+1-j
          if(scala.ge.0.0) then
            write(6,5) jwr,(nint(a(i,jwr)/scale),i=ib,ie,iskp)
    5       format(1x,i3,24i5)
          else
            write(6,6) jwr,(a(i,jwr),i=ib,ie,iskp)
    6       format(1x,i2,12(e10.2))
          endif
        end do
C
        write(6,7)
    7   format(//)
C
      end do
C
      return
C
      end
C
      subroutine prxyz(label,time,a,im,iskp,jm,jskp,kb,ko,nko,scala)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes horizontal layers of a 3-D field with        *
C *                integers or floating point numbers.                 *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jskp ........ skipping interval for j               *
C *                ko .......... 1-D array of k-indices for output     *
C *                nko ......... number of elements in ko              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                                                                    *
C *                (NOTE that this combines functions of old prxyz and *
C *                 eprxyz)                                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm,kb
      real a(im,jm,kb)
      real time,scala
      integer ko(*)
      integer iskp,jskp,nko
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,jwr,k,iko,cols
C
      if(scala.ge.0.0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.0) scale = 1.0
      if (scala.eq.0.0) then
        amx=1.e-12
        do iko=1,nko
          k=ko(iko)
          do j=1,jm,jskp
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.0**(int(log10(amx)+100.0)-103)
        endif
      if(scala.gt.0.0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do iko=1,nko
C
        k=ko(iko)
C
        write(6,3) k
    3   format(3x,/' Layer k = ',i2)
C
        do ib=1,im,cols*iskp
C
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
C
          if(scala.ge.0.0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,2x,24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,12i10,/)
          endif
C
          do j=1,jm,jskp
            jwr=jm+1-j
            if(scala.ge.0.0) then
              write(6,6) jwr,(nint(a(i,jwr,k)/scale),i=ib,ie,iskp)
    6         format(1x,i3,24i5)
            else
              write(6,7) jwr,(a(i,jwr,k),i=ib,ie,iskp)
    7         format(1x,i2,12(e10.2))
            endif
          end do
C
          write(6,8)
    8     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine prxz(label,time,a,im,iskp,jm,kb,jo,njo,scala,dt,zz)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
C *                x- or i-direction .                                 *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                iskp ........ skipping interval for i               *
C *                jo .......... 1-D array of j-indices for output     *
C *                njo ......... number of elements in jo              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                dt(im,jm) ... total depth                           *
C *                zz(kb) ...... sigma coordinate                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer im,jm,kb
      real a(im,jm,kb),dt(im,jm),zz(kb)
      real time,scala
      integer jo(*)
      integer iskp,njo
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,k,ijo,cols
C
      if(scala.ge.0.0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.0) scale = 1.0
      if (scala.eq.0.0) then
        amx=1.e-12
        do  k=1,kb
          do ijo=1,njo
            j=jo(ijo)
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.0**(int(log10(amx)+100.0)-103)
        endif
      if(scala.gt.0.0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do ijo=1,njo
C
        j=jo(ijo)
C
        write(6,3) j
    3   format(3x,/' Section j =',i3)
C
        do ib=1,im,cols*iskp
C
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
C
          if(scala.ge.0.0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,'    i =  ',24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,'    i =  ',12i10,/)
          endif
C
          write(6,6) (nint(dt(i,j)),i=ib,ie,iskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
C
          do k=1,kb
            if(scala.ge.0.0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),i=ib,ie,iskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),i=ib,ie,iskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
C
          write(6,9)
    9     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine pryz(label,time,a,im,jm,jskp,kb,io,nio,scala,dt,zz)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
C *                y- or j-direction.                                  *
C *                                                                    *
C *                label ....... label for output                      *
C *                time ........ time (days)                           *
C *                a(im,jm,kb).. array to be printed                   *
C *                jskp ........ skipping interval for j               *
C *                io .......... 1-D array of i-indices for output     *
C *                nio ......... number of elements in io              *
C *                scala ....... < 0 for floating point numbers output *
C *                              0 for integer output, divisor for a   *
C *                                based on magnitudes of |a| values   *
C *                              > 0 for integer output, divisor for a *
C *                                given by scala                      *
C *                dt(im,jm) ... total depth                           *
C *                zz(kb) ...... sigma coordinate                      *
C *                                                                    *
C **********************************************************************
C
      implicit none
      integer im,jm,kb
      real a(im,jm,kb),dt(im,jm),zz(kb)
      real time,scala
      integer io(*)
      integer jskp,nio
      character label*(*)
      real amx,scale
      integer i,j,jb,je,k,iio,cols
C
      if(scala.ge.0.0) then
        cols=24
      else
        cols=12
      endif
C
      if (scala.lt.0.0) scale = 1.0
      if (scala.eq.0.0) then
        amx=1.e-12
        do  k=1,kb
          do j=1,jm,jskp
            do iio=1,nio
              i=io(iio)
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.0**(int(log10(amx)+100.0)-103)
        endif
      if(scala.gt.0.0) scale=scala
C
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
C
      do iio=1,nio
C
        i=io(iio)
C
        write(6,3) i
    3   format(3x,/' Section i =',i3)
C
        do jb=1,jm,cols*jskp
C
          je=jb+(cols-1)*jskp
          if(je.gt.jm) je=jm
C
          if(scala.ge.0.0) then
            write(6,4) (j,j=jb,je,jskp)
    4       format(/,'    j =  ',24i5,/)
          else
            write(6,5) (j,j=jb,je,jskp)
    5       format(/,'    j =  ',12i10,/)
          endif
C
          write(6,6) (nint(dt(i,j)),j=jb,je,jskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
C
          do k=1,kb
            if(scala.ge.0.0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),j=jb,je,jskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),j=jb,je,jskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
C
          write(6,9)
    9     format(//)
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine seamount
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up for seamount problem.                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real delh,delx,elejmid,elwjmid,ra,vel
      integer i,j,k
C
C     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
C
      delh=0.90
C
C     Grid size:
C
      delx=8000.0
C
C     Radius island or seamount:
C
      ra=25000.0
C
C     Current velocity:
C
      vel=0.20
C
C     Set up grid dimensions, areas of free surface cells, and
C     Coriolis parameter:
C
      do j=1,jm
        do i=1,im
C
C     For constant grid size:
C
C         dx(i,j)=delx
C         dy(i,j)=delx
C
C     For variable grid size:
C
          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.0
          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.0
C
          cor(i,j)=1.e-4
C
        end do
      end do
C
C     Calculate horizontal coordinates of grid points and rotation
C     angle.
C
C     NOTE that this is introduced solely for the benefit of any post-
C     processing software, and in order to conform with the requirements
C     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
C
C     There are four horizontal coordinate systems, denoted by the
C     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
C     "e" is an elevation point and "c" is a cell corner), as shown
C     below. In addition, "east_*" is an easting and "north_*" is a
C     northing. Hence the coordinates of the "u" points are given by
C     (east_u,north_u).
C
C     Also, if the centre point of the cell shown below is at
C     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
C     the coordinates of the western of the two "u" points,
C     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
C     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
C     coordinates of the southwestern corner point of the cell. The
C     southwest corner of the entire grid is at
C     (east_c(1,1),north_c(1,1)).
C
C                      |              |
C                    --c------v-------c--
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                      u      e       u
C                      |              |
C                      |              |
C                      |              |
C                      |              |
C                    --c------v-------c--
C                      |              |
C
C
C     NOTE that the following calculation of east_c and north_c only
C     works properly for a rectangular grid with east and north aligned
C     with i and j, respectively:
C
      do j=1,jm
        east_c(1,j)=0.0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
C
      do i=1,im
        north_c(i,1)=0.0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
C
C     The following works properly for any grid:
C
C     Elevation points:
C
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)
     $                  +east_c(i,j+1)+east_c(i+1,j+1))/4.0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)
     $                   +north_c(i,j+1)+north_c(i+1,j+1))/4.0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im-1
        east_e(i,jm)
     $    =((east_c(i,jm)+east_c(i+1,jm))*3.0
     $       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.0
        north_e(i,jm)
     $    =((north_c(i,jm)+north_c(i+1,jm))*3.0
     $       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.0
      end do
C
      do j=1,jm-1
        east_e(im,j)
     $    =((east_c(im,j)+east_c(im,j+1))*3.0
     $       -east_c(im-1,j)-east_c(im-1,j+1))/4.0
        north_e(im,j)
     $    =((north_c(im,j)+north_c(im,j+1))*3.0
     $       -north_c(im-1,j)-north_c(im-1,j+1))/4.0
      end do
C
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)
     $               -(east_e(im-2,jm)+east_e(im,jm-2))/2.0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)
     $               -(north_e(im-2,jm)+north_e(im,jm-2))/2.0
C
C     u-points:
C
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.0
        end do
      end do
C
C     Extrapolate ends:
C
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.0-east_c(i,jm-1))/2.0
        north_u(i,jm)=(north_c(i,jm)*3.0-north_c(i,jm-1))/2.0
      end do
C
C     v-points:
C
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.0
        end do
      end do
C
C     Extrapolate ends:
C
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.0-east_c(im-1,j))/2.0
        north_v(im,j)=(north_c(im,j)*3.0-north_c(im-1,j))/2.0
      end do
C
C     rot is the angle (radians, anticlockwise) of the i-axis relative
C     to east, averaged to a cell centre:
C
C     (NOTE that the following calculation of rot only works properly
C     for this particular rectangular grid)
C
      do j=1,jm
        do i=1,im
          rot(i,j)=0.0
        end do
      end do
C
C     Define depth:
C
      do i=1,im
        do j=1,jm
C
          h(i,j)=4500.0*(1.0-delh
     $                          *exp(-((east_c(i,j)
     $                                   -east_c((im+1)/2,j))**2
     $                                +(north_c(i,j)
     $                                   -north_c(i,(jm+1)/2))**2)
     $                                /ra**2))
          if(h(i,j).lt.1.0) h(i,j)=1.0
C
        end do
      end do
C
C     Close the north and south boundaries to form a channel:
C
      do i=1,im
        h(i,1)=1.0
        h(i,jm)=1.0
      end do
C
C     Calculate areas and masks:
C
      call areas_masks
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
      if(slmax.lt.1.0) call slpmax
C
C     Set initial conditions:
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=5.0+15.0*exp(zz(k)*h(i,j)/1000.0)-tbias
            sb(i,j,k)=35.0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
            ub(i,j,k)=vel*dum(i,j)
          end do
        end do
      end do
C
C     Initialise uab and vab as necessary
C     (NOTE that these have already been initialised to zero in the
C     main program):
C
      do j=1,jm
        do i=1,im
          uab(i,j)=vel*dum(i,j)
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
      do j=1,jm
        do i=1,im
C     No conditions necessary for this problem
        end do
      end do
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in the seamount problem, the east and west boundaries are open,
C     while the south and north boundaries are closed through the
C     specification of the masks fsm, dum and dvm):
C
      rfe=1.0
      rfw=1.0
      rfn=1.0
      rfs=1.0
C
      do j=2,jmm1
        uabw(j)=uab(2,j)
        uabe(j)=uab(imm1,j)
C
C     Set geostrophically conditioned elevations at the boundaries:
C
        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
      end do
C
C     Adjust boundary elevations so that they are zero in the middle
C     of the channel:
C
      elejmid=ele(jmm1/2)
      elwjmid=elw(jmm1/2)
      do j=2,jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j)
        elw(j)=(elw(j)-elwjmid)*fsm(2,j)
      end do
C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
C
      end
C
      subroutine slpmax
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Limits the maximum of:                              *
C *                                                                    *
C *                  <difference of depths>/<sum of depths>            *
C *                                                                    *
C *                for two adjacent cells. The maximum possible value  *
C *                is unity.                                           *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real mean,del
      integer i,j,loop
C
      do loop=1,10
C
C     Sweep right:
C
        do j=2,jm-1
C
          do i=2,im-1
            if(fsm(i,j).ne.0.0.and.fsm(i+1,j).ne.0.0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.0+del)
                h(i,j)=mean*(1.0-del)
              endif
            endif
          end do
C
C    Sweep left:
C
          do i=im-1,2,-1
            if(fsm(i,j).ne.0.0.and.fsm(i+1,j).ne.0.0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.0+del)
                h(i,j)=mean*(1.0-del)
              endif
            endif
          end do
C
        end do
C
C   Sweep up:
C
        do i=2,im-1
C
          do j=2,jm-1
            if(fsm(i,j).ne.0.0.and.fsm(i,j+1).ne.0.0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.0+del)
                h(i,j)=mean*(1.0-del)
              endif
            endif
          end do
C
C   Sweep down:
C
          do j=jm-1,2,-1
            if(fsm(i,j).ne.0.0.and.fsm(i,j+1).ne.0.0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.0+del)
                h(i,j)=mean*(1.0-del)
              endif
            endif
          end do
C
        end do
C
      end do
C
      return
C
      end
C
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates the antidiffusive velocity used to       *
C *                reduce the numerical diffusion associated with the  *
C *                upstream differencing scheme.                       *
C *                                                                    *
C *                This is based on a subroutine of Gianmaria Sannino  *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                The calculations have been simplified by removing   *
C *                the shock switch option.                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real ff(im,jm,kb)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      real sw
      real mol,abs_1,abs_2
      real value_min,epsilon
      real udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
C
      parameter (value_min=1.e-9,epsilon=1.0e-14)
C
C     Apply temperature and salinity mask:
C
      do k=1,kb
        do i=1,im
          do j=1,jm
            ff(i,j,k)=ff(i,j,k)*fsm(i,j)
          end do
        end do
      end do
C
C     Recalculate mass fluxes with antidiffusion velocity:
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.0
            else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.0
     $              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))
     $             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.0
            end if
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.0
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.0
     $             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))
     $            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.0
            end if
          end do
        end do
      end do
C
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.0
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))
     $             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.0
            end if
          end do
        end do
      end do
C
      return
C
      end
C
      subroutine vertvl(xflux,yflux)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Calculates vertical velocity.                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k.c'
C
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
C
C     Reestablish boundary conditions:
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.250*(dy(i,j)+dy(i-1,j))
     $                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do
C
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.250*(dx(i,j)+dx(i,j-1))
     $                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do
C
C     NOTE that, if one wishes to include freshwater flux, the
C     surface velocity should be set to vflux(i,j). See also
C     change made to 2-D volume conservation equation which
C     calculates elf.
C
        do j=2,jmm1
          do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
          end do
        end do
C
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)
     $                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)
     $                        +yflux(i,j+1,k)-yflux(i,j,k))
     $                        /(dx(i,j)*dy(i,j))
     $                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do
C
      return
C
      end
C
C      include 'pom2k.n'                                       ! *netCDF*
C
C     End of source code
C
C-----------------------------------------------------------------------
C
