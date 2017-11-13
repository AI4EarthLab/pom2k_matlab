#define INCVARDEF include 'pom2008.c'
#define ZWADMODIFY_V real:: dummy##__line__=wadmodify
#define ZWADMODIFY_C wadmodify=__line__
      program pom2008
      
!{                                                                       !
! ---------------------------------------------------------------------!
!     For recent changes search for !lyo:_20080415:                    !
!                                                                      !
!     Main changes:                                                    !
!                                                                      !
!     (1). George Mellor found bugs in subr.profq (also in pom2k); use !
!          the new one in this version.                                !
!     (2). nsmolar=1 case (for wet-and-dry runs with nwad=1) is        !
!          incomplete, and should NOT be used                          !
!                                                                      !
!                       ... l.oey (Apr/15/2008)                        !
!}                                                                      !
! ---------------------------------------------------------------------!
!     This POM code has wetting and drying (WAD) capability            !
!                                                                      !
!     Set nwad=1 for WAD,set =0 to turn off WAD                        !
!     Set BOTH nwad=0 & nsmolar=0 to recover the original non-WAD POM  !
!                                                                      !
!     For details see:                                                 !
!                                                                      !
!     O2005 = Oey [2005; OceanModelling,  9, 133-150]                  !
!     O2006 = Oey [2006; OceanModelling, 13, 176-195]                  !
!     O2007 = Oey, Ezer et al [2007; Ocean Dynamics, 57, 205-221]      !
!                                                                      !
!     or, http://www.aos.princeton.edu/WWWPUBLIC/PROFS/                !
!     also, .../WWWPUBLIC/PROFS/pom98_wad_release_descriptions.html    !
!                                                                      !
!     For changes that I made, search for the following keyword:       !
!                                                                      !
      ! lyo:!wad:  or simply "wad:"                                     !
      ! tne:!wad:  (changes taken from T.Ezer's pom98_wad seamount case)!
!                                                                      !
!     Three (3) files are needed to run the model test cases:          !
!                                                                      !
!     pom2008.f (code), pom2008.c (common blocks etc) & pom2008.n (netcdf)   !
!                                                                      !
!     see the runscript runpom2008* provided with these files.           !
!                                                                      !
!     "DOUBLE  PRECISION (i.e. REAL*8)" option can be used to          !
!     compile, with or without WAD; e.g. w/Intel ifort, use -r8:       !
!                                                                      !
!     ifort pom2008.f -o a.out -r8                                       !
!                                                                      !
!     the -r8 option is the same as specifying -real_size 64 or        !
!     -auto_double.  However, the code works with or without -r8       !
!                                                                      !
!     A link to the netCDF library also needs to be specified. For this!
!     and other details, see the runscript:                            !
!                                                                      !
!     runpom2008*                                                        !
!                                                                      !
!     NOTE: for nsmolar=1, the PARAMETER (IM=???,JM=???) in subroutines!
!     wadadvt2d & wadsmoladif must match that specified in pom2008.c     !
!     This is done using "include 'grid'"                              !
!                                                                      !
!     Send bugs and/or improvements to: lyo@princeton.edu              !
!                                                                      !
!                       ... l.oey (May/02/2007)                        !
!                                 (Jul/24,30/2007)                     !
!                                 (Aug/12,18/2007)                     !
!                                                                      !
!     Support from the Minerals Management Service (MMS) for making    !
!     the development of WAD in POM possible is acknowledged.          !
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
!                                                                      !
! lyo:Notes:                                                            !
!     This version incorporates WAD into pom2k.f, and also corrects    !
!     the bug found by Charles Tang, Glen Carter and Alain Caya (and   !
!     corrected in the 2006-05-03 version of pom2k.f).                 !
!     Search for "clyo:"for all the changes I have made.               !
!                                                                      !
!                       ... lyo (Jun/14/2006)                          !
!                                                                      !
! **********************************************************************
! *                                                                    *
! *   The last code change as rcorded in pom2k.change was on           *
! *                                                                    *
! *                     2006-05-03                                     *
! *                                  (adding IC from file)             *
! *                                                                    *
! * FUNCTION    :  This is a version of the three dimensional, time    *
! *                dependent, primitive equation, ocean model          *
! *                developed by Alan Blumberg and George Mellor with   *
! *                subsequent contributions by Leo Oey, Steve Brenner  *
! *                and others. It is now called the Princeton Ocean    *
! *                Model. Two references are:                          *
! *                                                                    *
! *                Blumberg, A.F. and G.L. Mellor; Diagnostic and      *
! *                  prognostic numerical circulation studies of the   *
! *                  South Atlantic Bight, J. Geophys. Res. 88,        *
! *                  4579-4592, 1983.                                  *
! *                                                                    *
! *                Blumberg, A.F. and G.L. Mellor; A description of a  *
! *                  three-dimensional coastal ocean circulation model,*
! *                  Three-Dimensional Coastal Ocean Models, Coastal   *
! *                  and Estuarine Sciences, 4, N.S. Heaps, ed.,       *
! *                  American Geophysical Union, 1-16, 1987.           *
! *                                                                    *
! *                In subroutine profq the model makes use of the      *
! *                turbulence closure sub-model described in:          *
! *                                                                    *
! *                Mellor, G.L. and T. Yamada; Development of a        *
! *                  turbulence closure model for geophysical fluid    *
! *                  problems, Rev. Geophys. Space Phys., 20, No. 4,   *
! *                  851-875, 1982.                                    *
! *            (note recent profq that includes breaking waves)        *
! *                                                                    *
! *                A user's guide is available:                        *
! *                                                                    *
! *                Mellor, G.L.; User's guide for a three-dimensional, *
! *                  primitive equation, numerical ocean model.        *
! *                  Princeton University Report, 1998.                *
! *                                                                    *
! *                In October 2001, the source code underwent          *
! *                revision by John Hunter of the University of        *
! *                Tasmania. Major aspects of the revision were:       *
! *                                                                    *
! *                (1) The revision was based on pom98 updated to      *
! *                    12/9/2001.                                      *
! *                (2) Declaration of all variables.                   *
! *                (3) Rationalisation of the input of all constants.  *
! *                (4) Modifications to the "printer" output.          *
! *                (5) Output to a netCDF file.                        *
! *                (6) Inclusion of surface freshwater flux.           *
! *                (7) Inclusion of atmospheric pressure.              *
! *                (8) Inclusion of an additional problem to check (6) *
! *                    and (7), above.                                 *
! *                (9) Inclusion of option for Smolarkiewicz           *
! *                    advection scheme.                               *
! *                                                                    *
! *                This revised version is functionally almost         *
! *                equivalent to pom98. The output to device 6 from    *
! *                the "seamount" problem should be almost the same,   *
! *                any differences being due to minor format changes   *
! *                and improvements in rounding.                       *
! *                                                                    *
! *                This revision was helped by the following people:   *
! *                Tal Ezer, Peter Holloway, George Mellor, Rich       *
! *                Signell, Ian Webster, Brian Williams and Emma Young.*
! *                                                                    *
! **********************************************************************
! *                                                                    *
! *                                  GENERAL NOTES                     *
! *                                                                    *
! *                1. All units are S.I. (M.K.S.) unless otherwise     *
! *                   stated. NOTE that time is in days from the start *
! *                   of the run.                                      *
! *                                                                    *
! *                2. "b", <nothing> and "f" refers to backward,       *
! *                   central and forward time levels.                 *
! *                                                                    *
! *                3. NetCDF output may be used. In order to omit/use  *
! *                   netCDF, comment/uncomment all statements         *
! *                   carrying the comment "*netCDF*" at the end of    *
! *                   the line (or set netcdf_file='nonetcdf')         *
! *                                                                    *
! *                4. NetCDF is version 3. An attempt has been made to *
! *                   conform to the NetCDF Climate and Forecast (CF)  *
! *                   Metadata Conventions, but this may not yet be    *
! *                   complete (see:                                   *
! *                                                                    *
! *          http://www.cgd.ucar.edu/cms/eaton/cf-metadata/index.html) *
! *                                                                    *
! *                5. In order to use netCDF, the program should be    *
! *                   compiled with the appropriate library. For       *
! *                   example, if using g77, you may need to type:     *
! *                                                                    *
! *                     g77 -o pom2k pom2k.f /usr/lib/libnetcdf.a      *
! *                                                                    *
! *                   You should also have the "include" file of       *
! *                   netCDF subroutines (pom2k.n).                    *
! *                                                                    *
! *                6. In order to use netCDF, you may need to change   *
! *                   the name of the "include" file in the statement: *
! *                                                                    *
! *                     include '/usr/include/netcdf.inc'              *
! *                                                                    *
! *                   in subroutine write_netcdf                       *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! *                                SOFTWARE LICENSING                  *
! *                                                                    *
! *                This program is free software; you can redistribute *
! *                it and/or modify it under the terms of the GNU      *
! *                General Public License as published by the Free     *
! *                Software Foundation, either Version 2 of the        *
! *                license, or (at your option) any later version.     *
! *                                                                    *
! *                This program is distributed in the hope that it     *
! *                will be useful, but without any warranty; without   *
! *                even the implied warranty of merchantability or     *
! *                fitness for a particular purpose. See the GNU       *
! *                General Public License for more details.            *
! *                                                                    *
! *                A copy of the GNU General Public License is         *
! *                available at http://www.gnu.org/copyleft/gpl.html   *
! *                or by writing to the Free Software Foundation, Inc.,*
! *                59 Temple Place - Suite 330, Boston, MA 02111, USA. *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
! lyo:!wad:Add WAD variables in pom2008.c
      INCVARDEF
! 
!     New declarations plus ispi,isp2i:
! 
      real aam_init,atot
      real cbcmax,cbcmin,darea
      real days,dte2,dvol
      real eaver
      real horcon
      real ispi,isp2i
      real period,prtd1,prtd2
      real saver,smoth,sw,swtch
      real taver,time0
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

      ZWADMODIFY_V
      real wadsmoth  ! lyo:!wad:
      real cflmin    ! lyo:_20080415:
      integer nsmolar  ! lyo:!wad:
! 
! ***********************************************************************
! 
!     source should agree with source_c in pom2008.c and source_n in
!     pom2008.n.
! 
      source='pom2008  2008-04-18'
! 
      if(source.ne.source_c) then
        write(6,7)
    7   format(/'Incompatible versions of program_ and include files ',    &
     &          '..... program terminated; in pom2008.f'/)
        stop
      endif
! 
! ***********************************************************************
! 
      small=1.e-8           ! Small value
! 
      pi=atan(1.e0)*4.e0    ! PI
! 
! ***********************************************************************
! 
!     Input of filenames and constants:
! 
!     NOTE that the array sizes im, jm and kb should be set in pom2008.c 
!          or the 'grid' file created by runpom2008
! 
! -----------------------------------------------------------------------
! 
      title='Run 1                                   ' ! run's title
! 
! -----------------------------------------------------------------------
! 
      netcdf_file='pom2008.nc'  ! netCDF output file
!     netcdf_file='nonetcdf'      ! disable netCDF output
! 
! -----------------------------------------------------------------------
! 
!     Problem number:
! 
!     iproblem      problem      initialisation
!                    type          subroutine
! 
!         1        seamount       seamount
! 
!         2        conservation   box
!                  box
! 
!         3        IC from file   file2ic
! 
!         4n       WAD prob#n, n=1,2 or 3     !lyo:!wad:
! 
      ZWADMODIFY_C
      iproblem=41
! 
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
!     WAD Control Parameters:                                          !
!                                                                      !
! lyo:!wad:set nwad=1 for WAD run, =0 to recover non-WAD run with min   !
!     depths set at say 10m (defined in routine called depending on    !
!                            value of iproblem)                        !
!                                                                      !
      nwad=1
!                                                                      !
! lyo:!wad:set nsmolar=1 if water depth is to be solved by Smolarkiewicz!
!     See O2005 or O2006.                                              !
!     NOTE: nwad=1 w/nsmolar=1 may require finer grid AND even a larger!
!     hc, below, say hc=0.1 instead of default hc=0.05                 !
!                                                                      !
      nsmolar=0
!                                                                      !
! lyo:!wad:Define hhi0, later hhi will be set to it:                    !
!                                                                      !
!     hhi = HIghest water [above MSL] ever expected, this should       !
!           be the observed_highest + say 1 meter to make sure;        !
!           here, MSL = Mean Sea Level
!                                                                      !
!     Choose hhi0=20 (below) assuming that the highest tide or surge   !
!     is not expected to rise above 20m wrt MSL                        !
!                                                                      !
      hhi0=20.e0
!                                                                      !
! lyo:!wad:Define hc in meters:                                         !
!                                                                      !
!     hc  = Thinnest fluid layer (i.e. water depth) below which cell   !
!           is considered dry; default is 0.05m but can be changed in  !
!           "params" below; used only if nwad=1                        !
!                                                                      !
      hc=0.05e0
!                                                                      !
! lyo:!wad:Set wadsmoth.gt.0 to (cosmetically) smooth out isolated,     !
!     thin-fluid wet cells. Smoothing is done if the wet cell is thin, !
!     thickness <= hc*(1+wadsmoth); recommended non-zero wadsmoth is   !
!     ~0.05 but it should not be >~0.1                                 !
      wadsmoth=0.e0
! 
! -----------------------------------------------------------------------
! 
      zsh=.01e0          ! Bottom Log-layer shift (metres) !lyo:!wad:
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Tidal amplitude for iproblem=41 w/ or w/o WAD:               !
!                                                                      !
! lyo:_20080415:
      tidamp=0.e0
      if (iproblem .eq. 41) tidamp=9.e0
! -----!WAD END---------------------------------------------------------!
!                                                                      !
! 
! -----------------------------------------------------------------------
! 
!       mode                     description
! 
!        2        2-D calculation (bottom stress calculated in advave)
! 
!        3        3-D calculation (bottom stress calculated in profu,v)
! 
!        4        3-D calculation with t and s held fixed
! 
      mode=3
! 
! -----------------------------------------------------------------------
! 
!     Advection scheme:
! 
!      nadv     Advection scheme
! 
!        1       Centred scheme, as originally provide in POM
!        2       Smolarkiewicz iterative upstream scheme, based on
!                subroutines provided by Gianmaria Sannino and Vincenzo
!                Artale
! 
      nadv=1
! 
! -----------------------------------------------------------------------
! 
!     Constants for Smolarkiewicz iterative upstream scheme.
! 
!     Number of iterations. This should be in the range 1 - 4. 1 is
!     standard upstream differencing; 3 adds 50% CPU time to POM:
! 
      ZWADMODIFY_C
      nitera=3  ! =2   !lyo:!wad:
! 
!     Smoothing parameter. This should preferably be 1, but 0 < sw < 1
!     gives smoother solutions with less overshoot when nitera > 1:
! 
      sw=1.e0  ! =0.5e0   !lyo:!wad:
! 
! -----------------------------------------------------------------------
! 
!     Index to indicate whether run to start from restart file
!     (nread=0: no restart input file; nread=1: restart input file):
! 
      nread=0
! 
! -----------------------------------------------------------------------
! 
!     External (2-D) time step (secs.) according to CFL:
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
! lyo:!wad:For 3-d runs w/rivers (strong stratification)                !
!     strong currents can (and physically they should) develop         !
!     in thin film of water over nearly dry cells, so dte and          !
!     isplit (below) may need to be smaller.  For dx~dy~1km,           !
!     I have used values as small as DTE=3.0 and ISPLIT=5.             !
!                                                                      !
!     Note#1: Since alpha=0 when nwad=1 (see below), the corresponding !
!             DTE needs to be smaller than the DTE used for nwad=0.    !
!     Note#2: ISPLIT should also be smaller than its value for no WAD; !
!             see Appendix A of Oey [2006; Ocean Modell. 13, 176-195]. !
!                                                                      !
!     For example, for the seamount problem that Tal Ezer has tested   !
!     (im,jm=65,49),  he found that DTE=60 & ISPLIT=30 worked for      !
!     nwad=0 but blew up for nwad=1. However, by decreasing ISPLIT     !
!     =5->10, the model with nwad=1 works for DTE=20->60, except that  !
!     isolated cells w/thin fluid layers (~5cm) appear, surrounded by  !
!     cells which are mostly dry; this probably is to be expected for  !
!     the relatively large dx&dy >2~4km used.  Using smaller DTE=15    !
!     & ISPLIT=5, or DTE=6 & ISPLIT=10 (as used below) can remove the  !
!     isolated thin-fluid cells.                                       !
! ---------------------------------------------------------------------!
! 
      dte=6.e0
! 
! -----------------------------------------------------------------------
! 
!     <Internal (3-D) time step>/<External (2-D) time step>
!     (dti/dte; dimensionless):
! 
      isplit=10 ! tne:!wad:dte=6.0 & isplit=30 work w/nwad=0 & hmax=200m
! 
! -----------------------------------------------------------------------
! 
!     Date and time of start of initial run of model in format (i.e.
!     UDUNITS convention)
! 
!       YYYY-MM-DD HH:MM:SS <+/->HH:MM
! 
!     where "<+/->HH:MM" is the time zone (positive eastwards from
!     Coordinated Universal Time). NOTE that the climatological time
!     axis (i.e. beginning of year zero, which does not exist in the
!     real-world calendar) has been used here. Insert your own date
!     and time as required:
! 
      time_start='2000-01-01 00:00:00 +00:00'
! 
! -----------------------------------------------------------------------
! 
      ZWADMODIFY_C
      days=0.025e0       ! run duration in days
! 
! -----------------------------------------------------------------------
! 
      ZWADMODIFY_C
      prtd1=0.0125e0     ! Initial print interval (days)
! 
! -----------------------------------------------------------------------
! 
      prtd2=1.e0         ! Final print interval (days)
! 
! -----------------------------------------------------------------------
! 
      swtch=1000.e0      ! Time to switch from prtd1 to prtd2
! 
! -----------------------------------------------------------------------
! 
      iskp=1             ! Printout skip interval in i
! 
! -----------------------------------------------------------------------
! 
      jskp=1             ! Printout skip interval in j
! 
! -----------------------------------------------------------------------
! 
!     Logical for inertial ramp (.true. if inertial ramp to be applied
!     to wind stress and baroclinic forcing, otherwise .false.)
! 
      lramp=.false.
! 
! -----------------------------------------------------------------------
! 
!     Reference density (recommended values: 1025 for seawater,
!     1000 for freswater; S.I. units):
! 
      rhoref=1025.e0
! 
! -----------------------------------------------------------------------
! 
      tbias=0.e0         ! Temperature bias (deg. C)
! 
! -----------------------------------------------------------------------
! 
      sbias=0.e0         ! Salinity bias
! 
! -----------------------------------------------------------------------
! 
      grav=9.806e0       ! gravity constant (S.I. units)
! 
! -----------------------------------------------------------------------
! 
      kappa=0.4e0        ! von Karman's constant
! 
! -----------------------------------------------------------------------
! 
      z0b=.01e0          ! Bottom roughness (metres)
! 
! -----------------------------------------------------------------------
! 
      cbcmin=.0025e0     ! Minimum bottom friction coeff.
! 
! -----------------------------------------------------------------------
! 
      cbcmax=1.e0        ! Maximum bottom friction coeff.
! 
! -----------------------------------------------------------------------
! 
! lyo:!wad:I would use horcon=0.1 (un-related to WAD)
      horcon=0.2e0       ! Smagorinsky diffusivity coeff.
! 
! -----------------------------------------------------------------------
! 
!     Inverse horizontal turbulent Prandtl number
!     (ah/am; dimensionless):
! 
!     NOTE that tprni=0.e0 yields zero horizontal diffusivity!
! 
      tprni=.2e0
! 
! -----------------------------------------------------------------------
! 
!     Background viscosity used in subroutines profq, proft, profu and
!     profv (S.I. units):
! 
      umol=2.e-5
! 
! -----------------------------------------------------------------------
! 
!     Maximum depth used in radiation boundary condition in subroutine
!     bcond (metres):
! 
      hmax=4500.e0
! 
! -----------------------------------------------------------------------
! 
!     Maximum magnitude of vaf (used in check that essentially tests
!     for CFL violation):
! 
      vmaxl=100.e0
! 
! -----------------------------------------------------------------------
! 
!     Maximum allowable value of:
! 
!       <difference of depths>/<sum of depths>
! 
!     for two adjacent cells (dimensionless). This is used in subroutine
!     slpmax. If >= 1, then slpmax is not applied:
! 
      slmax=2.e0
! 
! -----------------------------------------------------------------------
! 
!     Integers defining the number of logarithmic layers at the
!     surface and bottom (used by subroutine depth). The number of
!     logarithmic layers are kl1-2 at the surface and kb-kl2-1
!     at the bottom. For no log portions, set kl1=2 and kl2=kb-1:
! 
      kl1=6
      kl2=kb-2
! 
! -----------------------------------------------------------------------
! 
!     Water type, used in subroutine proft.
! 
!       ntp    Jerlov water type
! 
!        1            i
!        2            ia
!        3            ib
!        4            ii
!        5            iii
! 
      ntp=2
! 
! -----------------------------------------------------------------------
! 
!     Surface temperature boundary condition, used in subroutine proft:
! 
!       nbct   prescribed    prescribed   short wave
!              temperature      flux      penetration
! 
!        1        no           yes           no
!        2        no           yes           yes
!        3        yes          no            no
!        4        yes          no            yes
! 
      nbct=1
! 
! -----------------------------------------------------------------------
! 
!     Surface salinity boundary condition, used in subroutine proft:
! 
!       nbcs   prescribed    prescribed
!               salinity      flux
! 
!        1        no           yes
!        3        yes          no
! 
!     NOTE that only 1 and 3 are allowed for salinity.
! 
      nbcs=1
! 
! -----------------------------------------------------------------------
! 
!     Step interval during which external (2-D) mode advective terms are
!     not updated (dimensionless):
! 
      ispadv=5
! 
! -----------------------------------------------------------------------
! 
!     Constant in temporal filter used to prevent solution splitting
!     (dimensionless):
! 
      smoth=0.10e0
! 
! -----------------------------------------------------------------------
! 
!     Weight used for surface slope term in external (2-D) dynamic
!     equation (a value of alpha = 0.e0 is perfectly acceptable, but the
!     value, alpha=.225e0 permits a longer time step):
! 
      ZWADMODIFY_C
      alph0=0.225e0 ! lyo:!wad:use alph0 instead of alpha - defined later
! 
! -----------------------------------------------------------------------
! 
!     Initial value of aam:
! 
      aam_init=500.e0
! 
!     End of input of constants
! ***********************************************************************
! 
! --- Above are the default parameters, alternatively one can
! --- use parameters from a file created by runscript runpom2k_pow_wad !clyo:wad:
! lyo:wad:beg:
!     Overwrite some input constants: see "params" above in runpom2008
! lyo:wad:end:
      include 'params'
! 
! ***********************************************************************
! 
!     Calculate some constants:
! 
      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:                                                             !
!                                                                      !
!     Dry-cell T/S relaxation using leapfrog-trapezoidal:              !
      ZWADMODIFY_C
      dtrat=dti/(0.5*86400.)  ! 1/2 day relaxation time-scale
      cwetrlx1=(1.-dtrat)/(1.+dtrat); cwetrlx2=2.*dtrat/(1.+dtrat)
!                                                                      !
      alpha=alph0*float(1-nwad) ! lyo:!wad: must=0 for WAD
!                                                                      !
      hhi=hhi0*float(nwad)
!                                                                      !
      nsmolar=nsmolar*nwad
!                                                                      !
!     Check nsmolar:                  !lyo:_20080415:                  !
!                                                                      !
      if( (nsmolar.ne.0).and.(mode.ne.2) )then
      write(6,'('' Stopped; incorrect inputs of '')')
      write(6,'('' nsmolar   = '',i10)') nsmolar
      write(6,'('' mode      = '',i10)') mode
      write(6,'('' For nsmolar ne 0, the present version works  '')')
      write(6,'('' only if mode = 2. So either run in barotropic '')')
      write(6,'('' mode = 2, or set nsmolar=0 and resubmit '')')
      stop
      endif
!                                                                      !
!     Check iproblem & nwad compatibility:                             !
!                                                                      !
      if( (iproblem.ne.41).and.(nwad.eq.1) )then ! lyo:_20080415:
      write(6,'('' iproblem   = '',i10)') iproblem
      write(6,'('' nwad       = '',i10)') nwad
      write(6,'('' iproblem must = 41 for nwad = 1 '')') ! lyo:_20080415:
      write(6,'('' Stopped: iproblem & nwad are not compatible '')')
      stop
      endif
!                                                                      !
! -----!WAD END---------------------------------------------------------!
!                                                                      !
! 
      iend=max(nint(days*24.e0*3600.e0/dti),2)
      iprint=nint(prtd1*24.e0*3600.e0/dti)
      iswtch=nint(swtch*24.e0*3600.e0/dti)
! 
      ispi=1.e0/float(isplit)
      isp2i=1.e0/(2.e0*float(isplit))
! 
! -----------------------------------------------------------------------
! 
!     Print initial summary:
! 
      write(6,'(/,'' source   = '',a40)') source
      write(6,'('' title      = '',a40/)') title
      write(6,'('' iproblem   = '',i10)') iproblem
      write(6,'('' nwad       = '',i10)') nwad        ! lyo:!wad:
      write(6,'('' nsmolar    = '',i10)') nsmolar     ! lyo:!wad:
      write(6,'('' hhi        = '',f10.3)') hhi       ! lyo:!wad:
      write(6,'('' hc         = '',f10.3)') hc        ! lyo:!wad:
      write(6,'('' wadsmoth   = '',f10.3)') wadsmoth  ! lyo:!wad:
      write(6,'('' zsh        = '',f10.6)') zsh         ! lyo:!wad:
      write(6,'('' tidamp     = '',f10.4)') tidamp
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
! 
! -----------------------------------------------------------------------
! 
!     Initialise boundary arrays:
! 
      do i=1,im
        vabn(i)=0.e0
        vabs(i)=0.e0
        eln(i)=0.e0     ! lyo:!wad:note:keep=0.e0 inst. of -hhi
        els(i)=0.e0     ! .. redefined later in wadseamount etc.
        do k=1,kb
          vbn(i,k)=0.e0
          vbs(i,k)=0.e0
          tbn(i,k)=0.e0
          tbs(i,k)=0.e0
          sbn(i,k)=0.e0
          sbs(i,k)=0.e0
        end do
      end do
! 
      do j=1,jm
        uabe(j)=0.e0
        uabw(j)=0.e0
        ele(j)=0.e0     ! lyo:!wad:note:keep=0.e0 inst. of -hhi
        elw(j)=0.e0     ! .. redefined later in wadseamount etc.
        do k=1,kb
          ube(j,k)=0.e0
          ubw(j,k)=0.e0
          tbe(j,k)=0.e0
          tbw(j,k)=0.e0
          sbe(j,k)=0.e0
          sbw(j,k)=0.e0
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
!     Initialise 2-D and 3-D arrays for safety (this may be overwritten
!     later):
! 
      do j=1,jm
        do i=1,im
          uab(i,j)=0.e0
          vab(i,j)=0.e0
          elb(i,j)=0.e0     ! lyo:!wad:note:keep=0.e0 inst. of -hhi
          etb(i,j)=0.e0     ! .. redefined later in wadseamount etc.
          e_atmos(i,j)=0.e0 ! lyo:!wad:note:=0.e0 is correct
          vfluxb(i,j)=0.e0
          vfluxf(i,j)=0.e0
          wusurf(i,j)=0.e0
          wvsurf(i,j)=0.e0
          wtsurf(i,j)=0.e0
          wssurf(i,j)=0.e0
          swrad(i,j)=0.e0
          drx2d(i,j)=0.e0
          dry2d(i,j)=0.e0
        end do
      end do
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            ub(i,j,k)=0.e0
            vb(i,j,k)=0.e0
          end do
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
!     Set up sigma layers:
! 
      if(iproblem.ne.3) call depth
! 
! -----------------------------------------------------------------------
! 
!     Read in grid data, and initial and lateral boundary conditions:
! 
      if(iproblem.eq.1) then
        call seamount
      else if(iproblem.eq.2) then
        call box
      else if(iproblem.eq.3) then
        call file2ic
      else if(iproblem.eq.41) then  ! lyo:!wad:
      ZWADMODIFY_C
        call wadseamount
      else
        write(6,8)
    8   format(/' Invalid value of iproblem ..... program_ terminated'/)
        stop
      endif
! 
!     Inertial period for temporal filter:
! 
! lyo:_20080415:
      ZWADMODIFY_C
      if (cor(im/2,jm/2).ne.0.0) then
         period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0
      else
         period=1.0  ! set to 1day
      endif
! lyo:_20080415:
!     Make sure that wetmask = fsm for non-WAD runs:
      if (nwad.eq.0) then
         wetmask(:,:)=fsm(:,:)
      endif
! 
!     Initialise time:
! 
      time0=0.e0
      time=0.e0
! 
!     Initial conditions:
! 
!     NOTE that lateral thermodynamic boundary conditions are often set
!     equal to the initial conditions and are held constant thereafter.
!     Users can of course create variable boundary conditions.
! 
      do i=1,im
        do j=1,jm
          ua(i,j)=uab(i,j)
          va(i,j)=vab(i,j)
          el(i,j)=elb(i,j)  ! lyo:!wad:note:already defined
          et(i,j)=etb(i,j)  !          in wadseamount w/hhi etc
          etf(i,j)=et(i,j)
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          w(i,j,1)=vfluxf(i,j)
        end do
      end do
! 
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
      ZWADMODIFY_C
      if (horcon.le.0.0) aam(:,:,:)=-horcon  ! lyo:_20080415:
! 
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
! 
      call dens(s,t,rho)
! 
      call baropg
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do
! 
! lyo:!wad:cbc_change:begins:-------------------------------------------!
!                                                                      !
! ---------------------------------------------------------------------!
!     Calculate bottom friction coefficient:                           !
!     ... see O2006, Appendix A; also Mellor,JPO,2002:Osc.Blayer...    !
! ---------------------------------------------------------------------!
!                                                                      !
      ZWADMODIFY_C
      do j=1,jm
        do i=1,im
!         cbc(i,j)=(kappa/log((    (1.e0+zz(kbm1))*h(i,j))/z0b))**2
          cbc(i,j)=(kappa/log((zsh+(1.e0+zz(kbm1))*d(i,j))/z0b))**2
! lyo:!wad:lyo's originals:---------------------------------------------!
! vers.1:  cbc(i,j)=(kappa/log((zsh+(zz(kbm1)-z(kb))*d(i,j))/z0b))**2   !
! vers.2:  cbc(i,j)=max(cbcmin,                                         !
!    $        (kappa/log((zsh+(zz(kbm1)-z(kb))*d(i,j))/z0b))**2)       !
! lyo:!wad:lyo's originals:---------------------------------------------!
! 
          cbc(i,j)=max(cbcmin,cbc(i,j))
! 
!     If the following is invoked, then it is probable that the wrong
!     choice of z0b or vertical spacing has been made:
! 
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
! lyo:!wad:cbc_change:ends:---------------------------------------------!
! 
!     Calculate external (2-D) CFL time step:
! 
! lyo:_20080415:
      cflmin=1.e10
      ZWADMODIFY_C
      do j=1,jm
        do i=1,im
          tps(i,j)=0.5e0/sqrt(1.e0/dx(i,j)**2+1.e0/dy(i,j)**2)            &
     &               /sqrt(grav*(h(i,j)+small))*fsm(i,j)
          if (fsm(i,j).ne.0.0) cflmin=min(cflmin,tps(i,j)) ! lyo:_20080415:
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
!     The following data are needed for a seamless restart. if nread=1,
!     data had been created by a previous run (see write(71) at end of
!     this program). nread=0 denotes a first time run.
! 
      if(nread.eq.1)                                                      &
      ZWADMODIFY_C
     &  read(70) time0,                                                   &
     &           wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,       &
     &           utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,                       &
     &           adx2d,ady2d,advua,advva,                                 &
     &           km,kh,kq,l,q2,q2b,aam,q2l,q2lb                           &
     &          ,wetmask,wmarsh   ! lyo:!wad:
! 
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do
! 
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Update cbc:                                                  !
!                                                                      !
      ZWADMODIFY_C
      cbc(:,:)=cbcmin      ! lyo:!wad_20080415:
      if (mode.ne.2) then  ! lyo:!wad_20080415:
      do j=1,jm
        do i=1,im
!         cbc(i,j)=cbcmin  !lyo:_20080415:
          cbc(i,j)=(kappa/log((zsh+(1.e0+zz(kbm1))*d(i,j))/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
      endif
!                                                                      !
! ---------------------------------------------------------------------!
! 
      time=time0
! 
! -----------------------------------------------------------------------
! 
!     Print geometry and other initial fields (select statements as
!     desired):
! 
      call prxy('grid increment in x, dx                 ',               &
     &          time,dx ,im,iskp,jm,jskp,0.e0)
! 
      call prxy('grid increment in y, dy                 ',               &
     &          time,dy ,im,iskp,jm,jskp,0.e0)
! 
      call prxy('Easting of elevation points, east_e     ',               &
     &          time,east_e ,im,iskp,jm,jskp,0.e0)
! 
      call prxy('Northing of elevation points, north_e   ',               &
     &          time,north_e,im,iskp,jm,jskp,0.e0)
! 
      call prxy('Easting of cell corners, east_c         ',               &
     &          time,east_c ,im,iskp,jm,jskp,0.e0)
! 
      call prxy('Northing of cell corners, north_c       ',               &
     &          time,north_c,im,iskp,jm,jskp,0.e0)
! 
      call prxy('Rotation angle of x-axis wrt. east, rot ',               &
     &          time,rot,im,iskp,jm,jskp,0.e0)
! 
      call prxy('Undisturbed water depth, h              ',               &
     &          time,h  ,im,iskp,jm,jskp,1.e1)
! 
      call prxy('Land mask, fsm                          ',     ! lyo:!wad:&
     &          time,fsm,im,iskp,jm,jskp,1.e0)
! 
      call prxy('Wet & dry mask, wetmask                 ',     ! lyo:!wad:&
     &          time,wetmask,im,iskp,jm,jskp,1.e0)
! 
      call prxy('u-velocity mask, dum                    ',               &
     &          time,dum,im,iskp,jm,jskp,1.e0)
! 
      call prxy('v-velocity mask, dvm                    ',               &
     &          time,dvm,im,iskp,jm,jskp,1.e0)
! 
      write(6,'('' Min CFL Time-step, cflmin = '',f10.4)') cflmin ! lyo:_20080415:

      call prxy('External (2-D) CFL time step, tps       ',               &
     &          time,tps,im,iskp,jm,jskp,1.e0)
! 
      call prxy('Bottom friction coefficient, cbc        ',     ! lyo:!wad:&
     &          time,cbc,im,iskp,jm,jskp,0.e0)
! 
!     Set sections for output:
! 
      ko(1)=1
      ko(2)=kb/2
      ko(3)=kb-1
! 
      call prxyz('Horizontally-averaged rho, rmean        ',              &
     &           time,rmean,im,iskp,jm,jskp,kb,ko,3,1.e-5)
! 
!     Set sections for output:
! 
      jo(1)=1
      jo(2)=jm/2
      jo(3)=jm-1
! 
      call prxz('Horizontally-averaged rho, rmean        ',               &
     &          time,rmean,im,iskp,jm,kb,jo,3,1.e-5,dt,zz)
! 
!     Set sections for output:
! 
      io(1)=1
      io(2)=im/2
      io(3)=im-1
! 
      call pryz('Horizontally-averaged rho, rmean        ',               &
     &          time,rmean,im,jm,jskp,kb,io,3,1.e-5,dt,zz)
! 
! -----------------------------------------------------------------------
! 
!     Initial conditions:
! 
!     Select print statements in printall as desired:
! 
      call printall
! 
! -----------------------------------------------------------------------
! 
!     Initialise netCDF output and output initial set of data:
! 
        if(netcdf_file.ne.'nonetcdf') then
      call write_netcdf(netcdf_file,1)                        ! *netCDF*
      call write_netcdf(netcdf_file,2)                        ! *netCDF*
        endif
! 
! -----------------------------------------------------------------------
! 
      do 9000 iint=1,iend      !  Begin internal (3-D) mode
! 
        time=dti*float(iint)/86400.e0+time0
! 
        if(lramp) then
          ramp=time/period
          if(ramp.gt.1.e0) ramp=1.e0
        else
          ramp=1.e0
        endif
! 
!       write(6,2) mode,iint,time
!   2   format(' mode,iint,time =',2i5,f9.2)
! 
! -----------------------------------------------------------------------
! 
!     Set time dependent, surface and lateral boundary conditions.
!     The latter will be used in subroutine bcond. Users may
!     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
!     wssurf, swrad and vflux.
! 
!     Introduce simple wind stress. Value is negative for westerly or
!     southerly winds. The following wind stress has been tapered
!     along the boundary to suppress numerically induced oscilations
!     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
!     To make a healthy surface Ekman layer, it would be well to set
!     kl1=9.
! 
        do j=2,jmm1
          do i=2,imm1
! 
      if(iproblem.ne.3) then     ! constant wind read in file2ic
! 
      ZWADMODIFY_C
!           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1))
            wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))                   &
     &                    *.25e0*(dvm(i,j+1)+dvm(i-1,j+1)                 &
     &                          +dvm(i-1,j)+dvm(i,j))
! --- no wind ----
!           wusurf(i,j)=0.e0  !lyo:_20080415:
            wvsurf(i,j)=0.e0
       endif
            e_atmos(i,j)=0.e0
            vfluxf(i,j)=0.e0
! 
!     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
!     the sea surface. See calculation of elf(i,j) below and subroutines
!     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
!     is no net flow across lateral boundaries, the basin volume will be
!     constant; if also vflux(i,j).ne.0, then, for example, the average
!     salinity will change and, unrealistically, so will total salt.
! 
            w(i,j,1)=vfluxf(i,j)
! 
!     Set wtsurf to the sensible heat, the latent heat (which involves
!     only the evaporative component of vflux) and the long wave
!     radiation:
! 
            wtsurf(i,j)=0.e0
! 
!     Set swrad to the short wave radiation:
! 
            swrad(i,j)=0.e0
! 
!     To account for change in temperature of flow crossing the sea
!     surface (generally quite small compared to latent heat effect)
! 
            tatm=t(i,j,1)+tbias    ! an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias)
! 
!     Set the salinity of water vapor/precipitation which enters/leaves
!     the atmosphere (or e.g., an ice cover)
! 
            satm=0.e0
            wssurf(i,j)=            vfluxf(i,j)*(satm-s(i,j,1)-sbias)
! 
          end do
        end do
! 
! lyo:
!     call powdriver(iprint,nread,z0b,cbcmin,iend/iprint,fsm)
! 
! -----------------------------------------------------------------------
! 
!     Set lateral viscosity:
! 
!     If mode=2 then initial values of aam2d are used. If one wishes
!     to use Smagorinsky lateral viscosity and diffusion for an
!     external (2-D) mode calculation, then appropiate code can be
!     adapted from that below and installed just before the end of the
!     "if(mode.eq.2)" loop in subroutine advave.
! 
!     Calculate Smagorinsky lateral viscosity:
! 
!       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
!                                     +.5*(du/dy+dv/dx)**2) )
! 
        if(mode.ne.2) then
          call advct(a,c,ee)
          call baropg
! 
      ZWADMODIFY_C
          if (horcon.gt.0.0) then ! lyo:_20080415:
          do k=1,kbm1
            do j=2,jmm1
              do i=2,imm1
                aam(i,j,k)=horcon*dx(i,j)*dy(i,j)                         &
     &                      *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2     &
     &                            +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2     &
     &                      +.5e0*(.25e0*(u(i,j+1,k)+u(i+1,j+1,k)         &
     &                                   -u(i,j-1,k)-u(i+1,j-1,k))        &
     &                      /dy(i,j)                                      &
     &                      +.25e0*(v(i+1,j,k)+v(i+1,j+1,k)               &
     &                             -v(i-1,j,k)-v(i-1,j+1,k))              &
     &                      /dx(i,j)) **2)
              end do
            end do
          end do
          endif  ! if (horcon.gt.0.0) then !lyo:_20080415:
! 
!     Form vertical averages of 3-D fields for use in external (2-D)
!     mode:
! 
          do j=1,jm
            do i=1,im
              adx2d(i,j)=0.e0
              ady2d(i,j)=0.e0
              drx2d(i,j)=0.e0
              dry2d(i,j)=0.e0
              aam2d(i,j)=0.e0
            end do
          end do
! 
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
! 
          call advave(tps)
! 
          do j=1,jm
            do i=1,im
              adx2d(i,j)=adx2d(i,j)-advua(i,j)
              ady2d(i,j)=ady2d(i,j)-advva(i,j)
            end do
          end do
! 
        endif
! 
! lyo:beg:changes by Tang et al.
        do j=1,jm
          do i=1,im
            egf(i,j)=el(i,j)*ispi
          end do
        end do
! 
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
! lyo:end:
! 
! -----------------------------------------------------------------------
! 
! lyomoving:decide to skip wriv & wmarsh for now; wriv in particular
!     needs to make consistent with vflux


        do 8000 iext=1,isplit    ! Begin external (2-D) mode
! 
!         write(6,3) iext,time
!   3     format(' iext,time =',i5,f9.2)
! 
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Water depth can be solved by Smolarkiewicz to prevent oscil. !
!                                                                      !
!     nsmolar=1 case (for wet-and-dry runs with nwad=1) is             !
!     incomplete, and should NOT be used                               !
!                                                                      !
      ZWADMODIFY_C
      if (nsmolar.eq.1 .and. mode.eq.2) then  ! lyo:_20080415:
        dd0(:,:)=h(:,:)+ el(:,:)*fsm(:,:)
        ddb(:,:)=h(:,:)+elb(:,:)*fsm(:,:)
        tps(:,:)=0.0           ! aam2d(:,:)
        call wadadvt2d(ddb,dd0,ddf,nitera,sw,dx,dy,ua,va,art,aru,arv,    &
     &    fsm,dum,dvm,tps,dte2)
        elf(:,:)=ddf(:,:)-h(:,:)
      else ! lyo:!wad:!
          do j=2,jm
            do i=2,im
              fluxua(i,j)=.25e0*(d(i,j)+d(i-1,j))                         &
     &                     *(dy(i,j)+dy(i-1,j))*ua(i,j)
              fluxva(i,j)=.25e0*(d(i,j)+d(i,j-1))                         &
     &                     *(dx(i,j)+dx(i,j-1))*va(i,j)
            end do
          end do
! 
!     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
!     with pom98.f. See also modifications to subroutine vertvl.
! 
          do j=2,jmm1
            do i=2,imm1
              elf(i,j)=elb(i,j)                                           &
     &                  +dte2*(-(fluxua(i+1,j)-fluxua(i,j)                &
     &                          +fluxva(i,j+1)-fluxva(i,j))/art(i,j)      &
     &                          -vfluxf(i,j))
            end do
          end do
!                                                                      !
      endif ! lyo:!wad:
!                                                                      !
! ---------------------------------------------------------------------!
! 
          call bcond(1)
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Check wet/dry on ELF:                                        !
!                                                                      !
      ZWADMODIFY_C
      if (nwad.eq.1) then
        do j=1,jm; do i=1,im
          DWET(I,J)=H(I,J)+ELF(I,J)*fsm(i,j)
        enddo; enddo
        do j=1,jm; do i=1,im
          wetmask(i,j)=fsm(i,j)
          if (dwet(i,j).le.hc) wetmask(i,j)=0.0
        enddo; enddo
        if (wadsmoth.gt.0.0) then
!       Smooth isolated wet spots:                                       !
          do j=1,jm; do i=1,im
            if( (wetmask(i,j) .gt.                                            &
     &        (wetmask(i-1,j)+wetmask(i+1,j)+wetmask(i,j-1)+wetmask(i,j+1)))  &
     &        .and. (dwet(i,j).le.hc*(1.e0+wadsmoth)) ) then
              wetmask(i,j)=0.0
!             elf(i,j)=elb(i,j); dwet(i,j)=h(i,j)+elf(i,j)*fsm(i,j)
            endif
          enddo; enddo
        endif
      endif
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
          if(mod(iext,ispadv).eq.0) call advave(tps)
! 
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=adx2d(i,j)+advua(i,j)                              &
     &                  -aru(i,j)*.25e0                                   &
     &                    *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))           &
     &                     +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))  &
     &                  +.25e0*grav*(dy(i,j)+dy(i-1,j))                   &
     &                    *(d(i,j)+d(i-1,j))                              &
     &                    *((1.e0-2.e0*alpha)                             &
     &                       *(el(i,j)-el(i-1,j))                         &
     &                      +alpha*(elb(i,j)-elb(i-1,j)                   &
     &                             +elf(i,j)-elf(i-1,j))                  &
     &                      +e_atmos(i,j)-e_atmos(i-1,j))                 &
     &                  +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
            end do
          end do
! 
          do j=2,jmm1
            do i=2,im
              uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))             &
     &                    *aru(i,j)*uab(i,j)                              &
     &                  -4.e0*dte*uaf(i,j))                               &
     &                 /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))            &
     &                     *aru(i,j))
            end do
          end do
! 
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=ady2d(i,j)+advva(i,j)                              &
     &                  +arv(i,j)*.25e0                                   &
     &                    *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))           &
     &                     +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))  &
     &                  +.25e0*grav*(dx(i,j)+dx(i,j-1))                   &
     &                    *(d(i,j)+d(i,j-1))                              &
     &                    *((1.e0-2.e0*alpha)*(el(i,j)-el(i,j-1))         &
     &                      +alpha*(elb(i,j)-elb(i,j-1)                   &
     &                             +elf(i,j)-elf(i,j-1))                  &
     &                      +e_atmos(i,j)-e_atmos(i,j-1))                 &
     &                  +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
            end do
          end do
! 
          do j=2,jm
            do i=2,imm1
              vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))             &
     &                    *vab(i,j)*arv(i,j)                              &
     &                  -4.e0*dte*vaf(i,j))                               &
     &                 /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))            &
     &                     *arv(i,j))
            end do
          end do
! 
          call bcond(2)
! 
          if(iext.eq.(isplit-2))then
            do j=1,jm
              do i=1,im
                etf(i,j)=.25e0*smoth*elf(i,j)
              end do
            end do
! 
          else if(iext.eq.(isplit-1)) then
! 
            do j=1,jm
              do i=1,im
                etf(i,j)=etf(i,j)+.5e0*(1.-.5e0*smoth)*elf(i,j)
              end do
            end do
! 
          else if(iext.eq.isplit) then
! 
            do j=1,jm
              do i=1,im
                etf(i,j)=(etf(i,j)+.5e0*elf(i,j))*fsm(i,j)
              end do
            end do
! 
          endif
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Check wet/dry on UAF & VAF:                                  !
!                                                                      !
      ZWADMODIFY_C
      if (nwad.eq.1) then
        do j=1,jm; do i=2,im
          if (0.5*(dwet(i,j)+dwet(i-1,j)).le.hc) then
            uaf(i,j)=0.0
          else                      ! Set flux=0 if coming from dry cell:
            if (wetmask(i-1,j).eq.0.0 .and. uaf(i,j).gt.0.0) uaf(i,j)=0.0
            if (wetmask(i,j)  .eq.0.0 .and. uaf(i,j).lt.0.0) uaf(i,j)=0.0
          endif
        enddo; enddo
! 
        do j=2,jm; do i=1,im
          if (0.5*(dwet(i,j)+dwet(i,j-1)).le.hc) then
            vaf(i,j)=0.0
          else                      ! Set flux=0 if coming from dry cell:
            if (wetmask(i,j-1).eq.0.0 .and. vaf(i,j).gt.0.0) vaf(i,j)=0.0
            if (wetmask(i,j)  .eq.0.0 .and. vaf(i,j).lt.0.0) vaf(i,j)=0.0
          endif
        enddo; enddo
      endif
!                                                                      !
! ---------------------------------------------------------------------!
! 
!     Stop if velocity condition violated (generally due to CFL
!     criterion not being satisfied):
! 
          vamax=0.e0
! 
          do j=1,jm
            do i=1,im
              if(abs(vaf(i,j)).ge.vamax) then
                vamax=abs(vaf(i,j))
                imax=i
                jmax=j
              endif
            end do
          end do
! 
          if(vamax.le.vmaxl) then
! 
!     Apply filter to remove time split and reset time sequence:
! 
            do j=1,jm
              do i=1,im
                ua(i,j)=ua(i,j)                                           &
     &                   +.5e0*smoth*(uab(i,j)-2.e0*ua(i,j)+uaf(i,j))
                va(i,j)=va(i,j)                                           &
     &                   +.5e0*smoth*(vab(i,j)-2.e0*va(i,j)+vaf(i,j))
                el(i,j)=el(i,j)                                           &
     &                   +.5e0*smoth*(elb(i,j)-2.e0*el(i,j)+elf(i,j))
                elb(i,j)=el(i,j)
                el(i,j)=elf(i,j)
                d(i,j)=h(i,j)+el(i,j)
                uab(i,j)=ua(i,j)
                ua(i,j)=uaf(i,j)
                vab(i,j)=va(i,j)
                va(i,j)=vaf(i,j)
              end do
            end do
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Update cbc:                                                  !
!                                                                      !
      ZWADMODIFY_C
      if (mode.ne.2) then  ! lyo:_20080415:
      do j=1,jm
        do i=1,im
          cbc(i,j)=cbcmin
          cbc(i,j)=(kappa/log((zsh+(1.e0+zz(kbm1))*d(i,j))/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do
      endif
!                                                                      !
! ---------------------------------------------------------------------!
! 
! lyo:beg:changes by Tang et al.
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
clyo:end:
            endif
! 
          endif
! 
 8000   end do        ! End of external (2-D) mode
! 
! -----------------------------------------------------------------------
! 
        if(vamax.le.vmaxl) then
! 
!     Continue with internal (3-D) mode calculation:
! 
          if((iint.ne.1.or.time0.ne.0.e0).and.mode.ne.2) then
! 
!     Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
! 
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
                end do
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=2,im
                  u(i,j,k)=(u(i,j,k)-tps(i,j))+                           &
     &                     (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
                end do
              end do
            end do
! 
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
                end do
              end do
            end do
! 
            do k=1,kbm1
              do j=2,jm
                do i=1,im
                  v(i,j,k)=(v(i,j,k)-tps(i,j))+                           &
     &                     (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
                end do
              end do
            end do
! 
!     vertvl calculates w from u, v, dt (h+et), etf and etb:
! 
            call vertvl(a,c)
            call bcond(5)
! 
! 
            do k=1,kb
              do j=1,jm
                do i=1,im
                  uf(i,j,k)=0.e0
                  vf(i,j,k)=0.e0
                end do
              end do
            end do
! 
!     Calculate q2f and q2lf using uf, vf, a and c as temporary
!     variables:
! 
            call advq(q2b,q2,uf,a,c)
            call advq(q2lb,q2l,vf,a,c)
            call profq(a,c,tps,dtef)
            call bcond(6)
! 
            do k=1,kb
              do j=1,jm
                do i=1,im
                  q2(i,j,k)=q2(i,j,k)                                     &
     &                       +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)            &
     &                                    -2.e0*q2(i,j,k))
                  q2l(i,j,k)=q2l(i,j,k)                                   &
     &                       +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)           &
     &                                    -2.e0*q2l(i,j,k))
                  q2b(i,j,k)=q2(i,j,k)
                  q2(i,j,k)=uf(i,j,k)
                  q2lb(i,j,k)=q2l(i,j,k)
                  q2l(i,j,k)=vf(i,j,k)
                end do
              end do
            end do
! 
!     Calculate tf and sf using uf, vf, a and c as temporary variables:
! 
            if(mode.ne.4) then
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:The followings prepare Oey's [1996] and Oey & Chen's [1992]  !
!     implementations of river and temperature surface conditions      !
!     i.e. using wriv -- but not implemented yet; they are not         !
!     directly for WAD, but ssurfwet is required later for wad         !
!                                                                      !
      ZWADMODIFY_C
            do j=1,jm
               do i=1,im
                 ssurfwet(i,j)=SSURF(i,j)*(1.-WRIV(i,j)/( WRIV(i,j)+      &
     &                                                     1.E-28)  )
               enddo
            enddo
!                                                                      !
! 
              if(nadv.eq.1) then
! 
                call advt1(tb,t,tclim,uf,a,c)
                call advt1(sb,s,sclim,vf,a,c)
! 
              else if(nadv.eq.2) then
! 
                call advt2(tb,t,tclim,uf,a,c,nitera,sw)
                call advt2(sb,s,sclim,vf,a,c,nitera,sw)
! 
              else
! 
                write(6,9)
    9           format(/'Invalid value for nadv ..... ',                  &
     &                 'program terminated'/)
                stop
! 
              endif
! 
              call proft(uf,wtsurf,tsurf,nbct,tps)
              call proft(vf,wssurf,ssurf,nbcs,tps)
              call bcond(4)
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Prepare dry-cell T/S "initial conditions" for next time-step !
!                                                                      !
      ZWADMODIFY_C
           if (nwad.eq.1) then
              do k=1,kb-1; do j=1,jm; do i=1,im
                 if (wetmask(i,j).eq.0.0) then
                    uf(i,j,k)=(cwetrlx1*tb(i,j,k)+cwetrlx2*tsurf(i,j)     &
     &                        *fsm(i,j))
                    vf(i,j,k)=(cwetrlx1*sb(i,j,k)+cwetrlx2*ssurfwet(i,j)  &
     &                        *fsm(i,j))
                 endif
                 enddo; enddo; enddo
           endif
! ---------------------------------------------------------------------!
! 
              do k=1,kb
                do j=1,jm
                  do i=1,im
                    t(i,j,k)=t(i,j,k)                                     &
     &                        +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)            &
     &                                     -2.e0*t(i,j,k))
                    s(i,j,k)=s(i,j,k)                                     &
     &                        +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)            &
     &                                     -2.e0*s(i,j,k))
                    tb(i,j,k)=t(i,j,k)
                    t(i,j,k)=uf(i,j,k)
                    sb(i,j,k)=s(i,j,k)
                    s(i,j,k)=vf(i,j,k)
                  end do
                end do
              end do
! 
              call dens(s,t,rho)
! 
            endif
! 
!     Calculate uf and vf:
! 
            call advu
            call advv
            call profu
            call profv
            call bcond(3)
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Check wet/dry on UF & VF:                                    !
!                                                                      !
      ZWADMODIFY_C
           if (nwad.eq.1) then
              do k=1,kb-1
                 do j=1,jm; do i=2,im
                    if (wetmask(i,j)*wetmask(i-1,j).eq.0.0)uf(i,j,k)=0.0
                    enddo; enddo
                 do j=2,jm; do i=1,im
                    if (wetmask(i,j)*wetmask(i,j-1).eq.0.0)vf(i,j,k)=0.0
                    enddo; enddo
              enddo
           endif
! ----!WAD END----------------------------------------------------------!
! 
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)                                       &
     &                      +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k)
                end do
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  u(i,j,k)=u(i,j,k)                                       &
     &                      +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)              &
     &                                   -2.e0*u(i,j,k)-tps(i,j))
                end do
              end do
            end do
! 
            do j=1,jm
              do i=1,im
                tps(i,j)=0.e0
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  tps(i,j)=tps(i,j)                                       &
     &                      +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k)
                end do
              end do
            end do
! 
            do k=1,kbm1
              do j=1,jm
                do i=1,im
                  v(i,j,k)=v(i,j,k)                                       &
     &                      +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)              &
     &                                   -2.e0*v(i,j,k)-tps(i,j))
                end do
              end do
            end do
! 
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
! 
          endif
! 
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
! 
        endif
! 
! -----------------------------------------------------------------------
! 
!     Beginning of print section:
! 
        if(iint.ge.iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
! 
        if(mod(iint,iprint).eq.0.or.vamax.gt.vmaxl) then
! 
          write(6,4) time,iint,iext,iprint
    4     format(/                                                        &
     &    '**************************************************',           &
     &    '**************************************************',           &
     &    '*************************'//                                   &
     &    ' time =',f9.4,', iint =',i8,', iext =',i8,', iprint =',i8,//)
! 
!     Select print statements in printall as desired:
! 
          call printall
!      ZWADMODIFY_C
          call wadout  ! lyo:!wad:
! 
          vtot=0.e0
          atot=0.e0
          taver=0.e0
          saver=0.e0
          eaver=0.e0
          do k=1,kbm1
            do j=1,jm
              do i=1,im
      ZWADMODIFY_C
                darea=dx(i,j)*dy(i,j)*wetmask(i,j) ! lyo:!wad:
                dvol=darea*dt(i,j)*dz(k)
                vtot=vtot+dvol
                taver=taver+tb(i,j,k)*dvol
                saver=saver+sb(i,j,k)*dvol
              end do
            end do
          end do
! 
          do j=1,jm
            do i=1,im
      ZWADMODIFY_C
              darea=dx(i,j)*dy(i,j)*wetmask(i,j) ! lyo:!wad:
              atot=atot+darea
              eaver=eaver+et(i,j)*darea
            end do
          end do
! 
          taver=taver/vtot
          saver=saver/vtot
          eaver=eaver/atot
          tsalt=(saver+sbias)*vtot
! 
          write(6,5) vtot,atot,eaver,taver,saver,tsalt
    5     format('vtot = ',e16.7,'   atot = ',e16.7,                      &
     &           '  eaver =',e16.7/'taver =',e16.7,                       &
     &           '   saver =',e16.7,'  tsalt =',e16.7)
! 
!     Write netCDF output:
! 
            if(netcdf_file.ne.'nonetcdf') then
          call write_netcdf(netcdf_file,2)                    ! *netCDF*
            endif
! 
          if(vamax.gt.vmaxl) then
! 
            write(6,4) time,iint,iext,iprint
! 
            call printall
! 
            write(6,6) vamax,imax,jmax
    6       format(///////////////////                                    &
     &             '************************************************'/    &
     &             '************ abnormal job end ******************'/    &
     &             '************* user terminated ******************'/    &
     &             '************************************************'/    &
     &             ' vamax =',e12.3,'   imax,jmax =',2i5)
! 
!     Close netCDF file:
! 
              if(netcdf_file.ne.'nonetcdf') then
            call write_netcdf(netcdf_file,3)                  ! *netCDF*
              endif
! 
            stop
! 
          endif
! 
        endif
! 
!     End of print section
! 
! -----------------------------------------------------------------------
! 
 9000 end do       !  End of internal (3-D) mode
! 
! -----------------------------------------------------------------------
! 
      write(6,4) time,iint,iext,iprint
! 
            call printall ! lyo:_20080415:final printing
! 
! 
!     Set levels for output:
! 
      ko(1)=1
      ko(2)=2
      ko(3)=kb/2
      ko(4)=kb-1
      ko(5)=kb
! 
      call prxyz('Vertical velocity, w                    ',              &
     &           time,w       ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
! 
!     call prxyz('Turbulent kinetic energy x 2, q2        ',
!    $           time,q2      ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
! 
!     Save this data for a seamless restart:
! 
! lyo:
!     call powsave
! 
      ZWADMODIFY_C
      write(71) time,                                                     &
     &  wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,                &
     &  utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,        &
     &  km,kh,kq,l,q2,q2b,aam,q2l,q2lb                                    &
     & ,wetmask,wmarsh   ! lyo:!wad:
! 
!     Close netCDF file:
! 
        if(netcdf_file.ne.'nonetcdf') then
      call write_netcdf(netcdf_file,3)                        ! *netCDF*
        endif
! 
! lyo:!wad:print final successful completion:
      write(6,10) time
   10 format(/2x,'JOB SUCCESSFULLY COMPLT.; time = ',1P1e13.5,' days')
! 
      stop
! 
      end program
! 
!     End of main program
! 
! -----------------------------------------------------------------------
! 
      subroutine advave(curv2d)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates horizontal advection and diffusion.      *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real curv2d(im,jm)
      integer i,j
! 
!     u-advection and diffusion:
! 
!     Advective fluxes:
! 
      do j=1,jm
        do i=1,im
          advua(i,j)=0.e0
        end do
      end do
! 
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.125e0*((d(i+1,j)+d(i,j))*ua(i+1,j)                 &
     &                       +(d(i,j)+d(i-1,j))*ua(i,j))                  &
     &                      *(ua(i+1,j)+ua(i,j))
        end do
      end do
! 
      do j=2,jm
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j)+d(i,j-1))*va(i,j)                   &
     &                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))            &
     &                      *(ua(i,j)+ua(i,j-1))
        end do
      end do
! 
!     Add viscous fluxes:
! 
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)                                         &
     &                 -d(i,j)*2.e0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))      &
     &                   /dx(i,j)
        end do
      end do
! 
      do j=2,jm
        do i=2,im
          tps(i,j)=.25e0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))            &
     &              *(aam2d(i,j)+aam2d(i,j-1)                             &
     &                +aam2d(i-1,j)+aam2d(i-1,j-1))                       &
     &              *((uab(i,j)-uab(i,j-1))                               &
     &                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))         &
     &               +(vab(i,j)-vab(i-1,j))                               &
     &                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25e0                        &
     &                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do
! 
      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)                            &
     &                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do
! 
!     u-advection and diffusion:
! 
      do j=1,jm
        do i=1,im
          advva(i,j)=0.e0
        end do
      end do
! 
!     Advective fluxes:
! 
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.125e0*((d(i,j)+d(i-1,j))*ua(i,j)                   &
     &                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))            &
     &                      *(va(i-1,j)+va(i,j))
        end do
      end do
! 
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j+1)+d(i,j))*va(i,j+1)                 &
     &                       +(d(i,j)+d(i,j-1))*va(i,j))                  &
     &                      *(va(i,j+1)+va(i,j))
        end do
      end do
! 
!     Add viscous fluxes:
! 
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)                                         &
     &                 -d(i,j)*2.e0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))      &
     &                   /dy(i,j)
        end do
      end do
! 
      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25e0                        &
     &                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do
! 
      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)                            &
     &                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do
! 
      if(mode.eq.2) then
! 
        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.5e0*(cbc(i,j)+cbc(i-1,j))                       &
     &                  *sqrt(uab(i,j)**2                                 &
     &                        +(.25e0*(vab(i,j)+vab(i,j+1)                &
     &                                 +vab(i-1,j)+vab(i-1,j+1)))**2)     &
     &                  *uab(i,j)
          end do
        end do
! 
        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.5e0*(cbc(i,j)+cbc(i,j-1))                       &
     &                  *sqrt(vab(i,j)**2                                 &
     &                        +(.25e0*(uab(i,j)+uab(i+1,j)                &
     &                                +uab(i,j-1)+uab(i+1,j-1)))**2)      &
     &                  *vab(i,j)
          end do
        end do
! 
        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.25e0                                             &
     &                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))      &
     &                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))     &
     &                   /(dx(i,j)*dy(i,j))
          end do
        end do
! 
        do j=2,jmm1
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25e0                          &
     &                  *(curv2d(i,j)*d(i,j)                              &
     &                    *(va(i,j+1)+va(i,j))                            &
     &                    +curv2d(i-1,j)*d(i-1,j)                         &
     &                    *(va(i-1,j+1)+va(i-1,j)))
          end do
        end do
! 
        do j=3,jmm1
          do i=2,imm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25e0                          &
     &                  *(curv2d(i,j)*d(i,j)                              &
     &                    *(ua(i+1,j)+ua(i,j))                            &
     &                    +curv2d(i,j-1)*d(i,j-1)                         &
     &                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
        end do
! 
      endif
! 
      return
! 
      end subroutine
! 
      subroutine advct(xflux,yflux,curv)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates the horizontal portions of momentum      *
! *                advection well in advance of their use in advu and  *
! *                advv so that their vertical integrals (created in   *
! *                the main program) may be used in the external (2-D) *
! *                mode calculation.                                   *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real curv(im,jm,kb)
      real dtaam
      integer i,j,k
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            curv(i,j,k)=0.e0
            advx(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.25e0*((v(i,j+1,k)+v(i,j,k))                      &
     &                         *(dy(i+1,j)-dy(i-1,j))                     &
     &                         -(u(i+1,j,k)+u(i,j,k))                     &
     &                         *(dx(i,j+1)-dx(i,j-1)))                    &
     &                       /(dx(i,j)*dy(i,j))
          end do
        end do
      end do
! 
!     Calculate x-component of velocity advection:
! 
!     Calculate horizontal advective fluxes:
! 
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.125e0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)           &
     &                           +(dt(i,j)+dt(i-1,j))*u(i,j,k))           &
     &                         *(u(i+1,j,k)+u(i,j,k))
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k)             &
     &                           +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))     &
     &                         *(u(i,j,k)+u(i,j-1,k))
          end do
        end do
      end do
! 
!    Add horizontal diffusive fluxes:
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k)=xflux(i,j,k)                                     &
     &                    -dt(i,j)*aam(i,j,k)*2.e0                        &
     &                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))         &
     &             *(aam(i,j,k)+aam(i-1,j,k)                              &
     &               +aam(i,j-1,k)+aam(i-1,j-1,k))
            yflux(i,j,k)=yflux(i,j,k)                                     &
     &                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))                 &
     &                            /(dy(i,j)+dy(i-1,j)                     &
     &                              +dy(i,j-1)+dy(i-1,j-1))               &
     &                            +(vb(i,j,k)-vb(i-1,j,k))                &
     &                            /(dx(i,j)+dx(i-1,j)                     &
     &                              +dx(i,j-1)+dx(i-1,j-1)))
! 
            xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)                         &
     &                          +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
          end do
        end do
      end do
! 
!     Do horizontal advection:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)                       &
     &                   +yflux(i,j+1,k)-yflux(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=3,imm1
            advx(i,j,k)=advx(i,j,k)                                       &
     &                   -aru(i,j)*.25e0                                  &
     &                     *(curv(i,j,k)*dt(i,j)                          &
     &                        *(v(i,j+1,k)+v(i,j,k))                      &
     &                       +curv(i-1,j,k)*dt(i-1,j)                     &
     &                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            advy(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do
! 
!     Calculate y-component of velocity advection:
! 
!     Calculate horizontal advective fluxes:
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*((dt(i,j)+dt(i-1,j))*u(i,j,k)             &
     &                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))     &
     &                         *(v(i,j,k)+v(i-1,j,k))
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k)=.125e0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)           &
     &                           +(dt(i,j)+dt(i,j-1))*v(i,j,k))           &
     &                         *(v(i,j+1,k)+v(i,j,k))
          end do
        end do
      end do
! 
!    Add horizontal diffusive fluxes:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))         &
     &             *(aam(i,j,k)+aam(i-1,j,k)                              &
     &               +aam(i,j-1,k)+aam(i-1,j-1,k))
            xflux(i,j,k)=xflux(i,j,k)                                     &
     &                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))                 &
     &                            /(dy(i,j)+dy(i-1,j)                     &
     &                              +dy(i,j-1)+dy(i-1,j-1))               &
     &                            +(vb(i,j,k)-vb(i-1,j,k))                &
     &                            /(dx(i,j)+dx(i-1,j)                     &
     &                              +dx(i,j-1)+dx(i-1,j-1)))
            yflux(i,j,k)=yflux(i,j,k)                                     &
     &                    -dt(i,j)*aam(i,j,k)*2.e0                        &
     &                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)
! 
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)                         &
     &                          +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do
! 
!     Do horizontal advection:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)                       &
     &                   +yflux(i,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=3,jmm1
          do i=2,imm1
            advy(i,j,k)=advy(i,j,k)                                       &
     &                   +arv(i,j)*.25e0                                  &
     &                     *(curv(i,j,k)*dt(i,j)                          &
     &                        *(u(i+1,j,k)+u(i,j,k))                      &
     &                       +curv(i,j-1,k)*dt(i,j-1)                     &
     &                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine advq(qb,q,qf,xflux,yflux)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates horizontal advection and diffusion, and  *
! *                vertical advection for turbulent quantities.        *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real qb(im,jm,kb),q(im,jm,kb),qf(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
! 
!     Do horizontal advection:
! 
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*(q(i,j,k)+q(i-1,j,k))                     &
     &                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.125e0*(q(i,j,k)+q(i,j-1,k))                     &
     &                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
! 
!     Do horizontal diffusion:
! 
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)                                     &
     &                    -.25e0*(aam(i,j,k)+aam(i-1,j,k)                 &
     &                            +aam(i,j,k-1)+aam(i-1,j,k-1))           &
     &                          *(h(i,j)+h(i-1,j))                        &
     &                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)         &
     &                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)                                     &
     &                    -.25e0*(aam(i,j,k)+aam(i,j-1,k)                 &
     &                            +aam(i,j,k-1)+aam(i,j-1,k-1))           &
     &                          *(h(i,j)+h(i,j-1))                        &
     &                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)         &
     &                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
! 
!     Do vertical advection, add flux terms, then step forward in time:
! 
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))       &
     &                 *art(i,j)/(dz(k)+dz(k-1))                          &
     &                 +xflux(i+1,j,k)-xflux(i,j,k)                       &
     &                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)                         &
     &                 *qb(i,j,k)-dti2*qf(i,j,k))                         &
     &                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine advt1(fb,f,fclim,ff,xflux,yflux)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Integrates conservative scalar equations.           *
! *                                                                    *
! *                This is centred scheme, as originally provide in    *
! *                POM (previously called advt).                       *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
! 
      do j=1,jm
        do i=1,im
           f(i,j,kb)=f(i,j,kbm1)
           fb(i,j,kb)=fb(i,j,kbm1)
        end do
      end do
! 
!     Do advective fluxes:
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*((dt(i,j)+dt(i-1,j))                       &
     &                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.25e0*((dt(i,j)+dt(i,j-1))                       &
     &                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do
! 
!     Add diffusive fluxes:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)                                     &
     &                    -.5e0*(aam(i,j,k)+aam(i-1,j,k))                 &
     &                         *(h(i,j)+h(i-1,j))*tprni                   &
     &                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)          &
     &                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)                                     &
     &                    -.5e0*(aam(i,j,k)+aam(i,j-1,k))                 &
     &                         *(h(i,j)+h(i,j-1))*tprni                   &
     &                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)          &
     &                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
! 
!     Do vertical advection:
! 
      do j=2,jmm1
        do i=2,imm1
          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
          zflux(i,j,kb)=0.e0
        end do
      end do
! 
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.5e0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do
! 
!     Add net horizontal fluxes and then step forward in time:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)                         &
     &                 +yflux(i,j+1,k)-yflux(i,j,k)                       &
     &                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
! 
            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)               &
     &                 -dti2*ff(i,j,k))                                   &
     &                 /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine advt2(fb,f,fclim,ff,xflux,yflux,nitera,sw)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Integrates conservative scalar equations.           *
! *                                                                    *
! *                This is a first-order upstream scheme, which        *
! *                reduces implicit diffusion using the Smolarkiewicz  *
! *                iterative upstream scheme with an antidiffusive     *
! *                velocity.                                           *
! *                                                                    *
! *                It is based on the subroutines of Gianmaria Sannino *
! *                (Inter-university Computing Consortium, Rome, Italy)*
! *                and Vincenzo Artale (Italian National Agency for    *
! *                New Technology and Environment, Rome, Italy),       *
! *                downloaded from the POM FTP site on 1 Nov. 2001.    *
! *                The calculations have been simplified by removing   *
! *                the shock switch option. It should be noted that    *
! *                this implementation does not include cross-terms    *
! *                which are in the original formulation.              *
! *                                                                    *
! *                fb,f,fclim,ff . as used in subroutine advt1         *
! *                xflux,yflux ... working arrays used to save memory  *
! *                nitera ........ number of iterations. This should   *
! *                                be in the range 1 - 4. 1 is         *
! *                                standard upstream differencing;     *
! *                                3 adds 50% CPU time to POM.         *
! *                sw ............ smoothing parameter. This should    *
! *                                preferably be 1, but 0 < sw < 1     *
! *                                gives smoother solutions with less  *
! *                                overshoot when nitera > 1.          *
! *                                                                    *
! *                Reference:                                          *
! *                                                                    *
! *                Smolarkiewicz, P.K.; A fully multidimensional       *
! *                  positive definite advection transport algorithm   *
! *                  with small implicit diffusion, Journal of         *
! *                  Computational Physics, 54, 325-362, 1984.         *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real sw
      integer nitera
      real fbmem(im,jm,kb),eta(im,jm)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      integer i,j,k,itera
! 
!     Calculate horizontal mass fluxes:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            xmassflux(i,j,k)=0.e0
            ymassflux(i,j,k)=0.e0
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j))                   &
     &                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do
! 
        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j))                   &
     &                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do
! 
      do j=1,jm
        do i=1,im
          fb(i,j,kb)=fb(i,j,kbm1)
          eta(i,j)=etb(i,j)
        end do
      end do
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            zwflux(i,j,k)=w(i,j,k)
            fbmem(i,j,k)=fb(i,j,k)
          end do
        end do
      end do
! 
!     Start Smolarkiewicz scheme:
! 
      do itera=1,nitera
! 
!     Upwind advection scheme:
! 
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.5e0                                          &
     &                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))    &
     &                        *fbmem(i-1,j,k)+                            &
     &                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))    &
     &                        *fbmem(i,j,k))
! 
              yflux(i,j,k)=0.5e0                                          &
     &                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))    &
     &                        *fbmem(i,j-1,k)+                            &
     &                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))    &
     &                        *fbmem(i,j,k))
            end do
          end do
        end do
! 
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,1)=0.e0
            if(itera.eq.1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
            zflux(i,j,kb)=0.e0
          end do
        end do
! 
        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.5e0                                          &
     &                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))          &
     &                       *fbmem(i,j,k)+                               &
     &                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))          &
     &                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do
! 
!     Add net advective fluxes and step forward in time:
! 
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)                       &
     &                 +yflux(i,j+1,k)-yflux(i,j,k)                       &
     &                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)          &
     &                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
! 
!     Calculate antidiffusion velocity:
! 
        call smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
! 
        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            do k=1,kb
              fbmem(i,j,k)=ff(i,j,k)
            end do
          end do
        end do
! 
!     End of Smolarkiewicz scheme
! 
      end do
! 
!     Add horizontal diffusive fluxes:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni         &
     &                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)                &
     &                   *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni         &
     &                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)                &
     &                   *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do
! 
!     Add net horizontal fluxes and step forward in time:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)         &
     &                               +yflux(i,j+1,k)-yflux(i,j,k))        &
     &                           /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine advu
! **********************************************************************
! *                                                                    *
! * ROUTINE NAME:  advu                                                *
! *                                                                    *
! * FUNCTION    :  Does horizontal and vertical advection of           *
! *                u-momentum, and includes coriolis, surface slope    *
! *                and baroclinic terms.                               *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer i,j,k
! 
!     Do vertical advection:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            uf(i,j,k)=0.e0
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.25e0*(w(i,j,k)+w(i-1,j,k))                         &
     &                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do
! 
!     Combine horizontal and vertical advection with coriolis, surface
!     slope and baroclinic terms:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)                                         &
     &                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)            &
     &                 -aru(i,j)*.25e0                                    &
     &                   *(cor(i,j)*dt(i,j)                               &
     &                      *(v(i,j+1,k)+v(i,j,k))                        &
     &                     +cor(i-1,j)*dt(i-1,j)                          &
     &                       *(v(i-1,j+1,k)+v(i-1,j,k)))                  &
     &                 +grav*.125e0*(dt(i,j)+dt(i-1,j))                   &
     &                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)        &
     &                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.e0)           &
     &                   *(dy(i,j)+dy(i-1,j))                             &
     &                 +drhox(i,j,k)
          end do
        end do
      end do
! 
!     Step forward in time:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))              &
     &                 *aru(i,j)*ub(i,j,k)                                &
     &                 -2.e0*dti2*uf(i,j,k))                              &
     &                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))             &
     &                  *aru(i,j))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine advv
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Does horizontal and vertical advection of           *
! *                v-momentum, and includes coriolis, surface slope    *
! *                and baroclinic terms.                               *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer i,j,k
! 
!     Do vertical advection:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            vf(i,j,k)=0.e0
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.25e0*(w(i,j,k)+w(i,j-1,k))                         &
     &                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do
! 
!     Combine horizontal and vertical advection with coriolis, surface
!     slope and baroclinic terms:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)                                         &
     &                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)            &
     &                 +arv(i,j)*.25e0                                    &
     &                   *(cor(i,j)*dt(i,j)                               &
     &                      *(u(i+1,j,k)+u(i,j,k))                        &
     &                     +cor(i,j-1)*dt(i,j-1)                          &
     &                       *(u(i+1,j-1,k)+u(i,j-1,k)))                  &
     &                 +grav*.125e0*(dt(i,j)+dt(i,j-1))                   &
     &                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)        &
     &                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.e0)           &
     &                   *(dx(i,j)+dx(i,j-1))                             &
     &                 +drhoy(i,j,k)
          end do
        end do
      end do
! 
!     Step forward in time:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))              &
     &                 *arv(i,j)*vb(i,j,k)                                &
     &                 -2.e0*dti2*vf(i,j,k))                              &
     &                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))             &
     &                  *arv(i,j))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine areas_masks
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates areas and masks.                         *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer i,j
! 
!     Calculate areas of "t" and "s" cells:
! 
      do j=1,jm
        do i=1,im
          art(i,j)=dx(i,j)*dy(i,j)
        end do
      end do
! 
!     Calculate areas of "u" and "v" cells:
! 
      do j=2,jm
        do i=2,im
          aru(i,j)=.25e0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25e0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
! 
      do j=1,jm
        aru(1,j)=aru(2,j)
        arv(1,j)=arv(2,j)
      end do
! 
      do i=1,im
        aru(i,1)=aru(i,2)
        arv(i,1)=arv(i,2)
      end do
! 
!     Initialise and set up free surface mask if no WAD: !lyo:!wad:
! 
      ZWADMODIFY_C
      if (nwad.eq.0) then ! lyo:!wad:
! 
!     do j=1,jm
!       do i=1,im
!         fsm(i,j)=0.e0
!         dum(i,j)=0.e0
!         dvm(i,j)=0.e0
!         if(h(i,j).gt.1.e0) fsm(i,j)=1.e0
!       end do
!     end do
! 
!     Set up velocity masks:
! 
!     do j=2,jm
!       do i=2,im
!         dum(i,j)=fsm(i,j)*fsm(i-1,j)
!         dvm(i,j)=fsm(i,j)*fsm(i,j-1)
!       end do
!     end do
! 
! lyo:_20080415:the above original pom2k version does NOT define dum & dvm
!     along i=j=1; use the folliwngs (from subroutine wadh).
      do j=1,jm; do i=1,im
        FSM(I,J)=1.; DUM(I,J)=1.; DVM(I,J)=1.
        IF(h(I,J).LE.1.0)THEN
          h(I,J)=1.0; FSM(I,J)=0.; DUM(I,J)=0.; DVM(I,J)=0.
        ENDIF
      enddo; enddo
      DO J=1,JM-1; DO I=1,IM
      IF(FSM(I,J).EQ.0..AND.FSM(I,J+1).NE.0.)DVM(I,J+1)=0.
      enddo; enddo
      DO J=1,JM;   DO I=1,IM-1
      IF(FSM(I,J).EQ.0..AND.FSM(I+1,J).NE.0.)DUM(I+1,J)=0.
      enddo; enddo
! 
      endif   ! if (nwad.eq.0) then       !lyo:!wad:
! 
      return
! 
      end subroutine
! 
      subroutine baropg
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates  baroclinic pressure gradient.           *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer i,j,k
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
          end do
        end do
      end do
! 
!     Calculate x-component of baroclinic pressure gradient:
! 
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))             &
     &                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do
! 
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)                                   &
     &                    +grav*.25e0*(zz(k-1)-zz(k))                     &
     &                      *(dt(i,j)+dt(i-1,j))                          &
     &                      *(rho(i,j,k)-rho(i-1,j,k)                     &
     &                        +rho(i,j,k-1)-rho(i-1,j,k-1))               &
     &                    +grav*.25e0*(zz(k-1)+zz(k))                     &
     &                      *(dt(i,j)-dt(i-1,j))                          &
     &                      *(rho(i,j,k)+rho(i-1,j,k)                     &
     &                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))                        &
     &                        *drhox(i,j,k)*dum(i,j)                      &
     &                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do
! 
!     Calculate y-component of baroclinic pressure gradient:
! 
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))             &
     &                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do
! 
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)                                   &
     &                    +grav*.25e0*(zz(k-1)-zz(k))                     &
     &                      *(dt(i,j)+dt(i,j-1))                          &
     &                      *(rho(i,j,k)-rho(i,j-1,k)                     &
     &                        +rho(i,j,k-1)-rho(i,j-1,k-1))               &
     &                    +grav*.25e0*(zz(k-1)+zz(k))                     &
     &                      *(dt(i,j)-dt(i,j-1))                          &
     &                      *(rho(i,j,k)+rho(i,j-1,k)                     &
     &                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))                        &
     &                        *drhoy(i,j,k)*dvm(i,j)                      &
     &                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do
! 
      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine bcond(idx)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Applies open boundary conditions.                   *
! *                                                                    *
! *                Closed boundary conditions are automatically        *
! *                enabled through specification of the masks, dum,    *
! *                dvm and fsm, in which case the open boundary        *
! *                conditions, included below, will be overwritten.    *
! *                                                                    *
! *                                The C-Grid                          *
! *                                **********                          *
! *                                                                    *
! *                The diagram below and some of the text was provided *
! *                by D.-S. Ko. It is for the case where u and v are   *
! *                the primary boundary conditions together with t and *
! *                s (co-located with e). e = el itself is rather      *
! *                unimportant and is substituted from an adjacent     *
! *                interior point.                                     *
! *                                                                    *
! *                The dotted line (.......) bounds the interior       *
! *                (non-boundary) grid points. In general, only those  *
! *                variables in the interior are computed and          *
! *                variables at open boundary have to be specified.    *
! *                All interpolations are centered in space except     *
! *                those at lateral open boundary where an upstream    *
! *                scheme is usually used. "BU" indictates a line of   *
! *                points where the boundary conditions should be      *
! *                specified.                                          *
! *                                                                    *
! *                Horizontal locations of e(el), t and s (etc.) are   *
! *                coincident. Unused points are indicated by "NU" and *
! *                "*". However, for attractive output, adjacent       *
! *                interior values may be filled in at these points.   *
! *                                                                    *
! *                People not acquainted with sigma coordinates have   *
! *                often asked what kind of boundary condition is      *
! *                applied along closed horizontal boundaries.         *
! *                Although the issue is not as important as it might  *
! *                be  for z-level grids, a direct answer is "half-    *
! *                slip" which, of course, is between free slip and    *
! *                non-slip.                                           *
! *                                                                    **********
! -------------------------------------------------------------------------------*
!     |                               N O R T H                                 *
!     |                                                                         *
!     |    1     2     3           i-1   i    i+1         im-2  im-1   im       *
! -----+-------------------------------------------------------------------------*
!     |   NU BC BC                                                    BC BC     *
!     |    v  v  v                                                     v  v     *
!     |                                                                         *
!     |BC> u* e  u  e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e  u  e  <BC*
!     |    |     |     |           |     |     |           |     |     |        *
!  jm |BC> +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v- <BC*
!     |    |     | ....|...........|.....|.....|...........|.....|.... |        *
!     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
!     |    |     | :   |           |     |     |           |     |   : |        *
! jm-1|    +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v-    *
!     |    |     | :   |           |     |     |           |     |   : |        *
!     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
!     |    |     | :   |           |     |     |           |     |   : |        *
! jm-2|    +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v-    *
!     |            :                                                 :          *
! W   |       .  . :.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .: .  .    E*
! E   |            :                    INTERIOR                     :         A*
! S   |       .    :.     .     .     .     .     .     .     .     .:    .    S*
! T   |            :                                                 :         T*
!     |    |     | :   |           |     |     |           |     |   : |        *
!     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
!     |    |     | :   |           |     |     |           |     |   : |        *
!   3 |    +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v-    *
!     |    |     | :   |           |     |     |           |     |   : |        *
!     |    u* e  u :e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e: u  e     *
!     |    |     | ....|...........|.....|.....|...........|.....|.... |        *
!   2 |BC> +--v--+--v--+--v--   .  +--v--+--v--+--v--   .  +--v--+--v--+--v- <BC*
!     |    |     |     |           |     |     |           |     |     |        *
!     |BC> u* e  u  e  u  e  .  .  u  e  u  e  u  e  .  .  u  e  u  e  u  e  <BC*
!     |    |     |     |           |     |     |           |     |     |      --*
!   1 |NU> +--v*-+--v*-+--v*-   .  +--v*-+--v*-+--v*-   .  +--v*-+--v*-+--v* <NU*
!     |                                                                         *
!     |    ^  ^  ^                                                     ^  ^     *
!     |   NU BC BC                                                    BC BC     *
! -----+-------------------------------------------------------------------------*
!     |    1     2     3           i-1   i     i+1         im-2  im-1  im       *
!     |                                                                         *
!     |                                S O U T H                                *
! -------------------------------------------------------------------------------*
!                                                                               *
! *******************************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer idx
      real ga,u1,wm
      integer i,j,k
      ZWADMODIFY_V
      real etide  ! lyo:!wad:add etide
! 
      if(idx.eq.1) then
! 
! -----------------------------------------------------------------------
! 
!     External (2-D) boundary conditions:
! 
!     In this example, the governing boundary conditions are a radiation
!     condition on uaf in the east and in the west, and vaf in the north
!     and south. The tangential velocities are set to zero on both
!     boundaries. These are only one set of possibilities and may not
!     represent a choice which yields the most physically realistic
!     result.
! 
!     Elevation (in this application, elevation is not a primary
!     boundary condition):
! 
        do j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
! 
        do i=1,im
          elf(i,1)=elf(i,2)
          elf(i,jm)=elf(i,jmm1)
        end do
! 
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
! 
        return
! 
      else if(idx.eq.2) then
! 
!     External (2-D) velocity:
! 
! tne:!wad:
! Very large (9m) tide-like BC (imposed by tidal inflow vel. not elev.)
!      note that ELW already include hhi (see SEAMOUNT)
! lyo:!wad:Specify time-dependent tidal tidal amplitude w/period=1day,
!         & amplitude=9m; note that tide is to be added to the geostro-
!         phically balanced ELW below.  So ETIDE is referenced wrt MSL:
!         [but perhaps it will be more direct to specify elevation??]
!        etide=-9.e0*sin(2.e0*pi*time/1.e0)
         etide=-tidamp*sin(2.e0*pi*time/1.e0) ! minus so ebb near time=0
! 
        do j=2,jmm1
! 
!     East:
! 
      ZWADMODIFY_C
! tne:!wad:- for wad replace H with D
          uaf(im,j)=uabe(j)                                               &
!    &               +rfe*sqrt(grav/h(imm1,j))                            &
     &               +rfe*sqrt(grav/d(imm1,j))                            &
     &                         *(el(imm1,j)-ele(j))
          uaf(im,j)=ramp*uaf(im,j)
          vaf(im,j)=0.e0
! 
!     West:
! 
! tne:!wad:- for wad replace H with D add tide
      ZWADMODIFY_C
          uaf(2,j)=uabw(j)                                                &
     &              -rfw*sqrt(grav/d(2,j))                                &
     &                        *(el(2,j)-elw(j)-etide)
          uaf(2,j)=ramp*uaf(2,j)
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0
! 
        end do
! 
        do i=2,imm1
! 
!     North:
! 
! tne:!wad:- for wad replace H with D
      ZWADMODIFY_C
          vaf(i,jm)=vabn(i)                                               &
     &               +rfn*sqrt(grav/d(i,jmm1))                            &
     &                         *(el(i,jmm1)-eln(i))
          vaf(i,jm)=ramp*vaf(i,jm)
          uaf(i,jm)=0.e0
! 
!     South:
! 
! tne:!wad:- for wad replace H with D
      ZWADMODIFY_C
          vaf(i,2)=vabs(i)                                                &
     &              -rfs*sqrt(grav/d(i,2))                                &
     &                        *(el(i,2)-els(i))
          vaf(i,2)=ramp*vaf(i,2)
          vaf(i,1)=vaf(i,2)
          uaf(i,1)=0.e0
! 
        end do
! 
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
! 
        return
! 
      else if(idx.eq.3) then
! 
! -----------------------------------------------------------------------
! 
!     Internal (3-D) boundary conditions:
! 
!     Velocity (radiation conditions; smoothing is used in the direction
!     tangential to the boundaries):
! 
        do k=1,kbm1
          do j=2,jmm1
! 
!     East:
! 
      ZWADMODIFY_C
!           ga=sqrt(h(im,j)/hmax)
            ga=sqrt(d(im,j)/hmax)  ! lyo:!wad:replace h w/d
            uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k)           &
     &                     +.25e0*u(imm1,j+1,k))                          &
     &                  +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k)      &
     &                    +.25e0*u(im,j+1,k))
            vf(im,j,k)=0.e0
! 
!     West:
! 
      ZWADMODIFY_C
            ga=sqrt(d(1,j)/hmax)  ! lyo:!wad:replace h w/d
            uf(2,j,k)=ga*(.25e0*u(3,j-1,k)+.5e0*u(3,j,k)                  &
     &                    +.25e0*u(3,j+1,k))                              &
     &                 +(1.e0-ga)*(.25e0*u(2,j-1,k)+.5e0*u(2,j,k)         &
     &                   +.25e0*u(2,j+1,k))
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0
          end do
        end do
! 
        do k=1,kbm1
          do i=2,imm1
! 
!     North:
! 
      ZWADMODIFY_C
            ga=sqrt(d(i,jm)/hmax)  ! lyo:!wad:replace h w/d
            vf(i,jm,k)=ga*(.25e0*v(i-1,jmm1,k)+.5e0*v(i,jmm1,k)           &
     &                     +.25e0*v(i+1,jmm1,k))                          &
     &                  +(1.e0-ga)*(.25e0*v(i-1,jm,k)+.5e0*v(i,jm,k)      &
     &                    +.25e0*v(i+1,jm,k))
            uf(i,jm,k)=0.e0
! 
!     South:
! 
      ZWADMODIFY_C
            ga=sqrt(d(i,1)/hmax)  ! lyo:!wad:replace h w/d
            vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k)                  &
     &                    +.25e0*v(i+1,3,k))                              &
     &                 +(1.e0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k)         &
     &                   +.25e0*v(i+1,2,k))
            vf(i,1,k)=vf(i,2,k)
            uf(i,1,k)=0.e0
          end do
        end do
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
! 
        return
! 
      else if(idx.eq.4) then
! 
!     Temperature and salinity boundary conditions (using uf and vf,
!     respectively):
! 
        do k=1,kbm1
          do j=1,jm
! 
!     East:
! 
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.e0) then
              uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
            else
              uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
              vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti                   &
     &              /((zz(k-1)-zz(k+1))*dt(imm1,j))
                uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
              endif
            endif
! 
!     West:
! 
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
            else
              uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
              vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti                         &
     &              /((zz(k-1)-zz(k+1))*dt(2,j))
                uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
              endif
            endif
          end do
        end do
! 
        do k=1,kbm1
          do i=1,im
! 
!     North:
! 
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
            else
              uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
              vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti                   &
     &              /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
              endif
            endif
! 
!     South:
! 
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
            else
              uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
              vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
              if(k.ne.1.and.k.ne.kbm1) then
                wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti                         &
     &              /((zz(k-1)-zz(k+1))*dt(i,2))
                uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
              endif
            endif
          end do
        end do
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
! 
        return
! 
      else if(idx.eq.5) then
! 
!     Vertical velocity boundary conditions:
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
! 
        return
! 
      else if(idx.eq.6) then
! 
!     q2 and q2l boundary conditions:
! 
        do k=1,kb
          do j=1,jm
! 
!     East:
! 
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
            if(u1.le.0.e0) then
              uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
            else
              uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
              vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
            endif
! 
!     West:
! 
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
            if(u1.ge.0.e0) then
              uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
            else
              uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
              vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
            endif
          end do
        end do
! 
        do k=1,kb
          do i=1,im
! 
!     North:
! 
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
            if(u1.le.0.e0) then
              uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
            else
              uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
              vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
            endif
! 
!     South:
! 
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
            if(u1.ge.0.e0) then
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
            else
              uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
              vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
            endif
          end do
        end do
! 
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do
! 
        return
! 
      endif
! 
      end subroutine
! 
      subroutine bcondorl(idx)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  This is an optional subroutine replacing  bcond and *
! *                using Orlanski's scheme (J. Comp. Phys. 21, 251-269,*
! *                1976), specialized for the seamount problem. To     *
! *                make it work for the seamount problem, I (G.M.)     *
! *                have had to add an extra condition on an "if"       *
! *                statement in the t and s open boundary conditions,  *
! *                which involves the sign of the normal velocity.     *
! *                Thus:                                               *
! *                                                                    *
! *            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k), *
! *                                                                    *
! *                plus 3 others of the same kind.                     *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer idx
      real cl,denom
      integer i,j,k
! 
      if(idx.eq.1) then
! 
! -----------------------------------------------------------------------
! 
!     External (2-D) boundary conditions:
! 
!     In this example the governing boundary conditions are a radiation
!     condition on uaf(im,j) in the east and an inflow uaf(2,j) in the
!     west. The tangential velocities are set to zero on both
!     boundaries. These are only one set of possibilities and may not
!     represent a choice which yields the most physically realistic
!     result.
! 
!     Elevation (in this application, elevation is not a primary
!     boundary condition):
! 
        do  j=1,jm
          elf(1,j)=elf(2,j)
          elf(im,j)=elf(imm1,j)
        end do
! 
        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do
! 
        return
! 
      else if(idx.eq.2) then
! 
!     External (2-D) velocity:
! 
        do j=2,jmm1
! 
!     West:
! 
          uaf(2,j)=ramp*uabw(j)-sqrt(grav/h(2,j))*(el(2,j)-elw(j))
          uaf(1,j)=uaf(2,j)
          vaf(1,j)=0.e0
! 
!     East:
! 
          uaf(im,j)=ramp*uabe(j)                                          &
     &               +sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
          vaf(im,j)=0.e0
! 
        end do
! 
        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do
! 
        return
! 
      else if(idx.eq.3) then
! 
! -----------------------------------------------------------------------
! 
!     Internal (3-D) boundary conditions:
! 
!     Eastern and western radiation boundary conditions according to
!     Orlanski's explicit scheme:
! 
        do k=1,kbm1
          do j=2,jmm1
! 
!     West:
! 
            denom=(uf(3,j,k)+ub(3,j,k)-2.e0*u(4,j,k))
            if(denom.eq.0.e0)denom=0.01e0
            cl=(ub(3,j,k)-uf(3,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(2,j,k)=(ub(2,j,k)*(1.e0-cl)+2.e0*cl*u(3,j,k))              &
     &                 /(1.e0+cl)
            uf(1,j,k)=uf(2,j,k)
            vf(1,j,k)=0.e0
! 
!     East:
! 
            denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.e0*u(im-2,j,k))
            if(denom.eq.0.e0)denom=0.01e0
            cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(im,j,k)=(ub(im,j,k)*(1.e0-cl)+2.e0*cl*u(im-1,j,k))         &
     &                  /(1.e0+cl)
            vf(im,j,k)=0.e0
          end do
        end do
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do
! 
        return
! 
      else if(idx.eq.4) then
! 
!     Temperature and salinity boundary conditions (using uf and vf,
!     respectively):
! 
        do k=1,kbm1
          do j=1,jm
! 
!     West:
! 
            ubw(j,k)=ub(2,j,k)
            denom=(uf(2,j,k)+tb(2,j,k)-2.e0*t(3,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(tb(2,j,k)-uf(2,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(1,j,k)=(tb(1,j,k)*(1.e0-cl)+2.e0*cl*t(2,j,k))/(1.e0+cl)
            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k)
! 
            denom=(vf(2,j,k)+sb(2,j,k)-2.e0*s(3,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(sb(2,j,k)-vf(2,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            vf(1,j,k)=(sb(1,j,k)*(1.e0-cl)+2.e0*cl*s(2,j,k))/(1.e0+cl)
            if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) vf(1,j,k)=sbw(j,k)
! 
!     East:
! 
            ube(j,k)=ub(im,j,k)
            denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.e0*t(im-2,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uf(im,j,k)=(tb(im,j,k)*(1.e0-cl)+2.e0*cl*t(im-1,j,k))         &
     &                  /(1.e0+cl)
            if(cl.eq.0.e0.and.ube(j,k).le.0.e0) uf(im,j,k)=tbe(j,k)
! 
            denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.e0*s(im-2,j,k))
            if(denom.eq.0.e0) denom=0.01e0
            cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            vf(im,j,k)=(sb(im,j,k)*(1.e0-cl)+2.e0*cl*s(im-1,j,k))         &
     &                  /(1.e0+cl)
            if(cl.eq.0.e0.and.ube(j,k).le.0.e0) vf(im,j,k)=sbe(j,k)
! 
          end do
        end do
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
! 
        return
! 
      else if(idx.eq.5) then
! 
!     Vertical velocity boundary conditions:
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do
! 
        return
! 
      else if(idx.eq.6) then
! 
!     q2 and q2l boundary conditions:
! 
        do k=1,kb
! 
          do j=1,jm
            uf(im,j,k)=1.e-10
            vf(im,j,k)=1.e-10
            uf(1,j,k)=1.e-10
            vf(1,j,k)=1.e-10
          end do
! 
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do
! 
        return
! 
      endif
! 
      end subroutine
! 
      subroutine box
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Sets up conservation box problem.                   *
! *                                                                    *
! *                This basin uses the same grid as the seamount       *
! *                problem, but it has a flat bottom, is surrounded by *
! *                walls and is initialised with uniform salinity and  *
! *                temperature. It is forced by a surface input of     *
! *                water of the same temperature and salinity as the   *
! *                water in the basin. Therefore, the temperature and  *
! *                salinity in the basin should not change, and the    *
! *                free surface should fall at a rate vflux. It is also*
! *                forced by a steady atmospheric pressure field which *
! *                depresses the southwestern half of the model by 1 m *
! *                and elevates the northeastern half of the model by  *
! *                1 m.                                                *
! *                                                                    *
! *                Since this problem defines its own fixed e_atmos,   *
! *                tatm, satm and e_atmos, comment out corresponding   *
! *                declarations after the do 9000 statement in main    *
! *                program.                                            *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real depth,delx,tatm,satm
      integer i,j,k
! 
!     Water depth:
! 
      depth=4500.e0
! 
!     Grid size:
! 
      delx=8000.e0
! 
!     Set up grid dimensions, areas of free surface cells, and
!     Coriolis parameter:
! 
      do j=1,jm
        do i=1,im
! 
!     For constant grid size:
! 
!         dx(i,j)=delx
!         dy(i,j)=delx
! 
!     For variable grid size:
! 
          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.e0
          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.e0
! 
          cor(i,j)=1.e-4
! 
        end do
      end do
! 
!     Calculate horizontal coordinates of grid points and rotation
!     angle.
! 
!     NOTE that this is introduced solely for the benefit of any post-
!     processing software, and in order to conform with the requirements
!     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
! 
!     There are four horizontal coordinate systems, denoted by the
!     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
!     "e" is an elevation point and "c" is a cell corner), as shown
!     below. In addition, "east_*" is an easting and "north_*" is a
!     northing. Hence the coordinates of the "u" points are given by
!     (east_u,north_u).
! 
!     Also, if the centre point of the cell shown below is at
!     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
!     the coordinates of the western of the two "u" points,
!     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
!     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
!     coordinates of the southwestern corner point of the cell. The
!     southwest corner of the entire grid is at
!     (east_c(1,1),north_c(1,1)).
! 
!                      |              |
!                    --c------v-------c--
!                      |              |
!                      |              |
!                      |              |
!                      |              |
!                      u      e       u
!                      |              |
!                      |              |
!                      |              |
!                      |              |
!                    --c------v-------c--
!                      |              |
! 
! 
!     NOTE that the following calculation of east_c and north_c only
!     works properly for a rectangular grid with east and north aligned
!     with i and j, respectively:
! 
      do j=1,jm
        east_c(1,j)=0.e0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
! 
      do i=1,im
        north_c(i,1)=0.e0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
! 
!     The following works properly for any grid:
! 
!     Elevation points:
! 
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)                          &
     &                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)                       &
     &                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im-1
        east_e(i,jm)                                                      &
     &    =((east_c(i,jm)+east_c(i+1,jm))*3.e0                            &
     &       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
        north_e(i,jm)                                                     &
     &    =((north_c(i,jm)+north_c(i+1,jm))*3.e0                          &
     &       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
      end do
! 
      do j=1,jm-1
        east_e(im,j)                                                      &
     &    =((east_c(im,j)+east_c(im,j+1))*3.e0                            &
     &       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
        north_e(im,j)                                                     &
     &    =((north_c(im,j)+north_c(im,j+1))*3.e0                          &
     &       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
      end do
! 
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)                       &
     &               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)                    &
     &               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
! 
!     u-points:
! 
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
! 
!     v-points:
! 
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
! 
!     rot is the angle (radians, anticlockwise) of the i-axis relative
!     to east, averaged to a cell centre:
! 
!     (NOTE that the following calculation of rot only works properly
!     for this particular rectangular grid)
! 
      do j=1,jm
        do i=1,im
          rot(i,j)=0.e0
        end do
      end do
! 
!     Define depth:
! 
      do i=1,im
        do j=1,jm
          h(i,j)=depth
        end do
      end do
! 
!     Close the north and south boundaries:
! 
      do i=1,im
        h(i,1)=1.e0
        h(i,jm)=1.e0
      end do
! 
!     Close the east and west boundaries:
! 
      do j=1,jm
        h(1,j)=1.e0
        h(im,j)=1.e0
      end do
! 
!     Calculate areas and masks:
! 
      call areas_masks
! 
!     Adjust bottom topography so that cell to cell variations
!     in h do not exceed parameter slmax:
! 
      if(slmax.lt.1.e0) call slpmax
! 
!     Set tbias and sbias here for test (tbias and sbias would
!     normally only be set in the main program):
! 
      tbias=10.e0
      sbias=20.e0
      write(6,1) tbias,sbias
    1 format(/' tbias and sbias changed in subroutine_ box to:'/           &
     &         2f10.3//)
! 
!     Set initial conditions:
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=20.e0-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
          end do
        end do
      end do
! 
!     Initialise uab and vab as necessary
!     (NOTE that these have already been initialised to zero in the
!     main program):
! 
      do j=1,jm
        do i=1,im
!     No conditions necessary for this problem
        end do
      end do
! 
!     Set surface boundary conditions, e_atmos, vflux, wusurf,
!     wvsurf, wtsurf, wssurf and swrad, as necessary
!     (NOTE:
!      1. These have all been initialised to zero in the main program.
!      2. The temperature and salinity of inflowing water must be
!         defined relative to tbias and sbias.):
! 
      do j=1,jm
        do i=1,im
           ZWADMODIFY_C
! 
! lyo:!wad:!pom2k_bug:tsurf and ssurf were never defined, but should be:
             tsurf(i,j)=tb(i,j,1)
             ssurf(i,j)=sb(i,j,1)
! 
          if(i+j-57.le.0) then
            e_atmos(i,j)=1.e0
          else
            e_atmos(i,j)=-1.e0
          endif
! 
!     Ensure atmospheric pressure cannot make water depth go negative:
! 
          e_atmos(i,j)=min(e_atmos(i,j),h(i,j))
! 
          vfluxf(i,j)=-0.0001e0
! 
!     See main program, just after "Begin numerical integration", for
!     an explanation of these terms:
! 
          tatm=20.e0
          satm=35.e0
! 
        end do
      end do
! 
!     Initialise elb, etb, dt and aam2d:
! 
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
! 
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      ZWADMODIFY_C
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
! ---------------------------------------------------------------------!
!                                                                      !
      call dens(sb,tb,rho)                                                                      ! 
! 
!     Generated horizontally averaged density field (in this
!     application, the initial condition for density is a function
!     of z (the vertical cartesian coordinate) -- when this is not
!     so, make sure that rmean has been area averaged BEFORE transfer
!     to sigma coordinates):
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
! 
!     Set lateral boundary conditions, for use in subroutine bcond
!     (in this problem, all lateral boundaries are closed through
!     the specification of the masks fsm, dum and dvm):
! 
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
! 
!     Set thermodynamic boundary conditions (for the seamount
!     problem, and other possible applications, lateral thermodynamic
!     boundary conditions are set equal to the initial conditions and
!     are held constant thereafter - users may, of course, create
!     variable boundary conditions):
! 
      do k=1,kbm1
! 
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
! 
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine dens(si,ti,rhoo)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates (density-1000.)/rhoref.                  *
! *                                                                    *
! *                (see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech.,  *
! *                609-611.)                                           *
! *                                                                    *
! *                ti is potential temperature                         *
! *                                                                    *
! *                If using 32 bit precision, it is recommended that   *
! *                cr,p,rhor,sr,tr,tr2,tr3 and tr4 be made double      *
! *                precision, and the "e"s in the constants be changed *
! *                to "d"s.                                            *
! *                                                                    *
! * NOTE: if pressure is not used in dens, buoyancy term (boygr)       *
! *       in profq must be changed (see note in profq)                 *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real si(im,jm,kb),ti(im,jm,kb),rhoo(im,jm,kb)
      real cr,p,rhor,sr,tr,tr2,tr3,tr4  ! ,hij !lyo:!wad:define hij
      integer i,j,k
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
! 
            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr
! 
!     Approximate pressure in units of bars:
! 
! 
! lyo:!wad:Keep "p" defined as in original pom w/o hhi (i.e. "h" defined
!         wrt MSL), but set p=0 for cells where input "h" indicates
!         "dry" (i.e. where h-hhi < 0); i.e. use "pdens":
! 
            ZWADMODIFY_C
            p=pdens(i,j,k)
!           p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5
! 
            rhor=-0.157406e0+6.793952e-2*tr                               &
     &            -9.095290e-3*tr2+1.001685e-4*tr3                        &
     &            -1.120083e-6*tr4+6.536332e-9*tr4*tr
! 
            rhor=rhor+(0.824493e0-4.0899e-3*tr                            &
     &            +7.6438e-5*tr2-8.2467e-7*tr3                            &
     &            +5.3875e-9*tr4)*sr                                      &
     &            +(-5.72466e-3+1.0227e-4*tr                              &
     &            -1.6546e-6*tr2)*abs(sr)**1.5                            &
     &            +4.8314e-4*sr*sr
! 
            cr=1449.1e0+.0821e0*p+4.55e0*tr-.045e0*tr2                    &
     &          +1.34e0*(sr-35.e0)
            rhor=rhor+1.e5*p/(cr*cr)*(1.e0-2.e0*p/(cr*cr))
! 
            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)
! 
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine depth
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Establishes the vertical sigma grid with log        *
! *                distributions at the top and bottom and a linear    *
! *                distribution in between. The number of layers of    *
! *                reduced thickness are kl1-2 at the surface and      *
! *                kb-kl2-1 at the bottom. kl1 and kl2 are defined in  *
! *                the main program. For no log portions, set kl1=2    *
! *                and kl2=kb-1.                                       *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real delz
      integer kdz(12)
      integer k
! 
      data kdz/1,1,2,4,8,16,32,64,128,256,512,1024/
! 
      z(1)=0.e0
! 
      do k=2,kl1
        z(k)=z(k-1)+kdz(k-1)
      end do
! 
      delz=z(kl1)-z(kl1-1)
! 
      do k=kl1+1,kl2
        z(k)=z(k-1)+delz
      end do
! 
      do k=kl2+1,kb
        dz(k)=float(kdz(kb-k+1))*delz/float(kdz(kb-kl2))
        z(k)=z(k-1)+dz(k)
      end do
! 
      do k=1,kb
        z(k)=-z(k)/z(kb)
      end do
! 
      do k=1,kb-1
        zz(k)=0.5e0*(z(k)+z(k+1))
      end do
! 
      zz(kb)=2.e0*zz(kb-1)-zz(kb-2)
! 
      do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
! 
      dz(kb)=0.e0
      dzz(kb)=0.e0
! 
      write(6,1)
    1 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
! 
      do k=1,kb
        write(6,2) k,z(k),zz(k),dz(k),dzz(k)
    2   format((' ',i5,4f10.3))
      end do
! 
      write(6,3)
    3 format(//)
! 
      return
! 
      end subroutine
! 
      subroutine findpsi
! **********************************************************************
! *                                                                    *
! * ROUTINE NAME:  findpsi                                             *
! *                                                                    *
! * FUNCTION    :  Calculates the stream function, first assuming      *
! *                zero on the southern boundary and then, using the   *
! *                values on the western boundary, the stream function *
! *                is calculated again. If the elevation field is near *
! *                steady state, the two calculations should agree;    *
! *                otherwise not.                                      *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      integer i,j
! 
      do j=1,jm
        do i=1,im
          psi(i,j)=0.e0
        end do
      end do
! 
!     Sweep northward:
! 
      do j=2,jmm1
        do i=2,im
          psi(i,j+1)=psi(i,j)                                             &
     &                +.25e0*uab(i,j)*(d(i,j)+d(i-1,j))                   &
     &                  *(dy(i,j)+dy(i-1,j))
        end do
      end do
! 
      call prxy('Streamfunction, psi from u              ',               &
     &          time,psi,im,iskp,jm,jskp,0.e0)
! 
!    Sweep eastward:
! 
      do j=2,jm
        do i=2,imm1
          psi(i+1,j)=psi(i,j)                                             &
     &                -.25e0*vab(i,j)*(d(i,j)+d(i,j-1))                   &
     &                  *(dx(i,j)+dx(i,j-1))
        end do
      end do
! 
      call prxy('Streamfunction, psi from v              ',               &
     &          time,psi,im,iskp,jm,jskp,0.e0)
! 
      return
! 
      end subroutine
! 
! ---------------------------------------------------------------------!
! lyo:_20080415:This version (of file2ic) is basically still the old one!
!     from pom2k - I have not converted it to setup for realistic      !
!     WAD runs - but see subr. wadseamount for an example of how-to.   !
! ---------------------------------------------------------------------!
      subroutine file2ic
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Sets up my own problem.                             *
! *                                                                    *
! * This example read IC from IC.dat file, generated by GRID.f in      *
! * GRID-DATA directory. Only minimal number of fields are read,       *
! * while others are calculated here.                                  *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real rad,re,dlat,dlon,cff
      integer i,j,k,m
      character*5 field
      rad=0.01745329
      re=6371.E3
! 
      write(6,'(/,'' Read grid and initial conditions '',/)')
! 
! --- 1D ---
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
! --- 2D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') east_e
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') north_e
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') h  ! lyo:!wad:note:read +- topography
! --- 3D ---
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') t
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') s
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') rmean
! --- Constant wind stress read here
! (for time dep. read in loop 9000 & interpolate in time)
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') wusurf
      read(40,'(a5)') field
      write(6,'(a5)') field
       read(40,'(8E12.5)') wvsurf
! 
! --- print vertical grid distribution
! 
      write(6,2)
    2 format(/2x,'k',7x,'z',9x,'zz',9x,'dz',9x,'dzz',/)
      write(6,'(''  '',/)')
      do k=1,kb
        write(6,3) k,z(k),zz(k),dz(k),dzz(k)
    3   format((' ',i5,4f10.3))
      end do
      write(6,'(''  '',//)')
! 
! lyo:_20080415:Move "calc. Coriolis etc" through "call areas_masks"
!     from below (after call dens) to here; otherwise, fsm will not
!     be defined but is used in subroutine dens. (This is a pom2k_bug).
! 
! --- calc. Coriolis Parameter
! 
        ZWADMODIFY_C
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29E-5*sin(north_e(i,j)*rad)
!           aam2d(i,j)=aam(i,j,1) !lyo:_20080415:pom2k_bug:aam not defined
            elb(i,j)=0.
            etb(i,j)=0.
            dt(i,j)=h(i,j)
          end do
        end do
! 
        do j=1,jm
          do i=2,im-1
            dx(i,j)=0.5*rad*re*sqrt(((east_e(i+1,j)-east_e(i-1,j))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i+1,j)-north_e(i-1,j))**2)
          end do
            dx(1,j)=dx(2,j)
            dx(im,j)=dx(im-1,j)
        end do
! 
        do i=1,im
          do j=2,jm-1
            dy(i,j)=0.5*rad*re*sqrt(((east_e(i,j+1)-east_e(i,j-1))
     1 *cos(north_e(i,j)*rad))**2+(north_e(i,j+1)-north_e(i,j-1))**2)
          end do
            dy(i,1)=dy(i,2)
            dy(i,jm)=dy(i,jm-1)
        end do
! 
!     Calculate areas and masks:
! 
      call areas_masks
! 
! --- calc. surface & lateral BC from climatology
! 
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
! 
!                    --- EAST & WEST BCs ---
        do j=1,jm
              ele(j)=0.
              elw(j)=0.
! --- other vel. BCs (fixed in time) can be specified here
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
!                    --- NORTH & SOUTH BCs ---
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
! 
!     Set initial conditions:
! 
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
! 
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      ZWADMODIFY_C
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
! ---------------------------------------------------------------------!
!                                                                      !
      call dens(sb,tb,rho)
      ZWADMODIFY_C ! deleted code here
! 
! --- the following grids are needed only for netcdf plotting
! 
!     Corner of cell points:
! 
      do j=2,jm
        do i=2,im
          east_c(i,j)=(east_e(i,j)+east_e(i-1,j)                          &
     &                  +east_e(i,j-1)+east_e(i-1,j-1))/4.e0
          north_c(i,j)=(north_e(i,j)+north_e(i-1,j)                       &
     &                   +north_e(i,j-1)+north_e(i-1,j-1))/4.e0
        end do
      end do
! 
! 
!     Extrapolate ends (approx.):
! 
      do i=2,im
        east_c(i,1)=2.*east_c(i,2)-east_c(i,3)
        north_c(i,1)=2.*north_c(i,2)-north_c(i,3)
      end do
        east_c(1,1)=2.*east_c(2,1)-east_c(3,1)
! 
      do j=2,jm
        east_c(1,j)=2.*east_c(2,j)-east_c(3,j)
        north_c(1,j)=2.*north_c(2,j)-north_c(3,j)
      end do
        north_c(1,1)=2.*north_c(1,2)-north_c(1,3)
! 
!     u-points:
! 
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
! 
!     v-points:
! 
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
! 
!     rot is the angle (radians, anticlockwise) of the i-axis relative
!     to east, averaged to a cell centre: (only needed for CDF plotting)
! 
      do j=1,jm
        do i=1,im-1
          rot(i,j)=0.
          dlat=north_e(i+1,j)-north_e(i,j)
          dlon= east_e(i+1,j)- east_e(i,j)
           if(dlon.ne.0.) rot(i,j)=atan(dlat/dlon)
        end do
       rot(im,j)=rot(im-1,j)
      end do
! 
!     Set lateral boundary conditions, for use in subroutine bcond
!     set all=0 for closed BCs.
!     Values=0 for vel BC only, =1 is combination of vel+elev.
      rfe=0.e0
      rfw=0.e0
      rfn=0.e0
      rfs=0.e0
! 
      return
      end subroutine
! 
      subroutine printall
! **********************************************************************
! *                                                                    *
! *                         POM2K SOURCE CODE                          *
! *                                                                    *
! * ROUTINE NAME:  printall                                            *
! *                                                                    *
! * FUNCTION    :  Prints a set of outputs to device 6                 *
! *                                                                    *
! *                Edit as approriate.                                 *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      integer io(100),jo(100),ko(100)
! 
      INCVARDEF
! 
!     2-D horizontal fields:
! 
          call prxy('Depth-averaged u, uab                   ',           &
     &              time,uab,im,iskp,jm,jskp,0.e0)
! 
          call prxy('Depth-averaged v, vab                   ',           &
     &              time,vab,im,iskp,jm,jskp,0.e0)
! 
          call prxy('Surface elevation, elb                  ',           &
     &              time,elb,im,iskp,jm,jskp,0.e0)
! 
!     Calculate and print streamfunction:
! 
          call findpsi
! 
          if(mode.ne.2) then
! 
!     2-D horizontal sections of 3-D fields:
! 
!     Set levels for output:
! 
            ko(1)=1
            ko(2)=kb/2
            ko(3)=kb-1
! 
            call prxyz('x-velocity, u                           ',        &
     &                 time,u    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
            call prxyz('y-velocity, v                           ',        &
     &                 time,v    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
            ko(1)=2
            call prxyz('z-velocity, w                           ',        &
     &                 time,w    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
            ko(1)=1
! 
            call prxyz('Potential temperature, t                ',        &
     &                 time,t    ,im,iskp,jm,jskp,kb,ko,3,1.e-2)
! 
            call prxyz('Salinity, s                              ',       &
     &                 time,s    ,im,iskp,jm,jskp,kb,ko,3,1.e-2)
! 
            call prxyz('(density-1000)/rhoref, rho              ',        &
     &                 time,rho  ,im,iskp,jm,jskp,kb,ko,3,1.e-5)
! 
!           call prxyz('Turbulent kinetic energy x 2, q2        ',
!    $                 time,q2   ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
!           call prxyz('Turbulent length scale, l               ',
!    $                 time,l    ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
            call prxyz('Horizontal kinematic viscosity, aam     ',        &
     &                 time,aam  ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
            call prxyz('Vertical kinematic viscosity, km        ',        &
     &                 time,km   ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
!           call prxyz('Vertical kinematic diffusivity, kh      ',
!    $                 time,kh   ,im,iskp,jm,jskp,kb,ko,3,0.e0 )
! 
!     Vertical sections of 3-D fields, normal to j-axis:
! 
!     Set sections for output:
! 
            jo(1)=1
            jo(2)=jm/2
            jo(3)=jm-1
! 
            call prxz('x-velocity, u                           ',         &
     &                time,u    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz)
! 
            call prxz('y-velocity, v                           ',         &
     &                time,v    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz)
! 
            call prxz('z-velocity, w                           ',         &
     &                time,w    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
! 
            call prxz('Potential temperature, t                ',         &
     &                time,t    ,im,iskp,jm,kb,jo,3,1.e-2,dt,zz)
! 
            call prxz('Salinity, s                             ',         &
     &                time,s    ,im,iskp,jm,kb,jo,3,1.e-2,dt,zz)
! 
            call prxz('(density-1000)/rhoref, rho              ',         &
     &                time,rho  ,im,iskp,jm,kb,jo,3,1.e-5,dt,zz)
! 
!           call prxz('Turbulent kinetic energy x 2, q2        ',
!    $                time,q2   ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
! 
!           call prxz('Turbulent length scale, l               ',
!    $                time,l    ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
! 
!           call prxz('Horizontal kinematic viscosity, aam     ',
!    $                time,aam  ,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz)
! 
!           call prxz('Vertical kinematic viscosity, km        ',
!    $                time,km   ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
! 
!           call prxz('Vertical kinematic diffusivity, kh      ',
!    $                time,kh   ,im,iskp,jm,kb,jo,3,0.e0 ,dt,z )
! 
!     Vertical sections of 3-D fields, normal to i-axis:
! 
!     Set sections for output:
! 
            io(1)=1
            io(2)=im/2
            io(3)=im-1
! 
            call pryz('x-velocity, u                           ',         &
     &                time,u    ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
! 
            call pryz('y-velocity, v                           ',         &
     &                time,v    ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
! 
            call pryz('z-velocity, w                           ',         &
     &                time,w    ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
! 
            call pryz('Potential temperature, t                ',         &
     &                time,t    ,im,jm,jskp,kb,io,3,1.e-2,dt,zz)
! 
!           call pryz('Salinity x rho / rhoref, s              ',
!    $                time,s    ,im,jm,jskp,kb,io,3,1.e-2,dt,zz)
! 
!           call pryz('(density-1000)/rhoref, rho              ',
!    $                time,rho  ,im,jm,jskp,kb,io,3,1.e-5,dt,zz)
! 
!           call pryz('Turbulent kinetic energy x 2, q2        ',
!    $                time,q2   ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
! 
!           call pryz('Turbulent length scale, l               ',
!    $                time,l    ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
! 
!           call pryz('Horizontal kinematic viscosity, aam     ',
!    $                time,aam  ,im,jm,jskp,kb,io,3,0.e0 ,dt,zz)
! 
!           call pryz('Vertical kinematic viscosity, km        ',
!    $                time,km   ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
! 
!           call pryz('Vertical kinematic diffusivity, kh      ',
!    $                time,kh   ,im,jm,jskp,kb,io,3,0.e0 ,dt,z )
! 
          endif
! 
      return
! 
      end subroutine
! 
      subroutine profq(sm,sh,dh,cc) ! lyo:_20080415:see notes at beg.
! **********************************************************************
! *                                        Updated: Sep. 24, 2003      *
! * FUNCTION    :  Solves for q2 (twice the turbulent kinetic energy), *
! *                q2l (q2 x turbulent length scale), km (vertical     *
! *                kinematic viscosity) and kh (vertical kinematic     *
! *                diffusivity), using a simplified version of the     *
! *                level 2 1/2 model of Mellor and Yamada (1982).      *
! * In this version, the Craig-Banner sub-model whereby breaking wave  *
! * tke is injected into the surface is included. However, we use an   *
! * analytical solution to the near surface tke equation to solve for  *
! * q2 at the surface giving the same result as C-B diffusion. The new *
! * scheme is simpler and more robust than the latter scheme.          *
! *                                                                    *
! * References                                                         *
! *   Craig, P. D. and M. L. Banner, Modeling wave-enhanced turbulence *
! *     in the ocean surface layer. J. Phys. Oceanogr., 24, 2546-2559, *
! *     1994.                                                          *
! *   Ezer, T., On the seasonal mixed-layer simulated by a basin-scale *
! *     ocean model and the Mellor-Yamada turbulence scheme,           *
! *     J. Geophys. Res., 105(C7), 16,843-16,855, 2000.                *
! *   Mellor, G.L. and T. Yamada, Development of a turbulence          *
! *     closure model for geophysical fluid fluid problems,            *
! *     Rev. Geophys. Space Phys., 20, 851-875, 1982.                  *
! *   Mellor, G. L., One-dimensional, ocean surface layer modeling,    *
! *     a problem and a solution. J. Phys. Oceanogr., 31(3), 790-809,  *
! *     2001.                                                          *
! *   Mellor, G.L. and A. Blumberg, Wave breaking and ocean surface    *
! *     thermal response, J. Phys. Oceanogr., 2003.                    *
! *   Stacey, M. W., Simulations of the wind-forced near-surface       *
! *     circulation in Knight Inlet: a parameterization of the         *
! *     roughness length. J. Phys. Oceanogr., 29, 1363-1367, 1999.     *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
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
! 
      integer i,j,k,ki
! 
      equivalence (prod,kn)
! 
      data a1,b1,a2,b2,c1/0.92e0,16.6e0,0.74e0,10.1e0,0.08e0/
      data e1/1.8e0/,e2/1.33e0/
      data sef/1.e0/
      data cbcnst/100./surfl/2.e5/shiw/0.0/
! 
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.e0*umol)*.5e0         &
     &                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.e0*umol)*.5e0         &
     &                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
!     The following section solves the equation:
! 
!       dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b
! 
!     Surface and bottom boundary conditions:
! 
      const1=(16.6e0**(2.e0/3.e0))*sef
! 
! initialize fields that are not calculated on all boundaries
! but are later used there
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
! 
      do j=1,jmm1
        do i=1,imm1
          utau2=sqrt((.5e0*(wusurf(i,j)+wusurf(i+1,j)))**2                &
     &                  +(.5e0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
! Wave breaking energy- a variant of Craig & Banner (1994)
! see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.e0
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2
! Surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2/grav
! 
          uf(i,j,kb)=sqrt((.5e0*(wubot(i,j)+wubot(i+1,j)))**2             &
     &                   +(.5e0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
! 
!    Calculate speed of sound squared:
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
! 
!     Calculate pressure in units of decibars:
! 
! 
! lyo:!wad:Keep "p" defined as in original pom w/o hhi (i.e. "h" defined
!         wrt MSL), but set p=0 for cells where input "h" indicates
!         "dry" (i.e. where h-hhi < 0); i.e. use "pdens":
! 
            ZWADMODIFY_C
            p=pdens(i,j,k)*10.e0
!           p=grav*rhoref*(-zz(k)* h(i,j))*1.e-4
            cc(i,j,k)=1449.1e0+.00821e0*p+4.55e0*tp -.045e0*tp**2         &
     &                 +1.34e0*(sp-35.0e0)
            cc(i,j,k)=cc(i,j,k)                                           &
     &                 /sqrt((1.e0-.01642e0*p/cc(i,j,k))                  &
     &                   *(1.e0-0.40e0*p/cc(i,j,k)**2))
          end do
        end do
      end do
! 
!     Calculate buoyancy gradient:
! 
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            ZWADMODIFY_C
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))                   &
! lyo:!wad:Note in MAIN, call profq is before call proft/dens. So T&S
!         (hence rho) are defined on sigma-grid using (h+etf),
!         so use "h+et" here:
!    $                    /(dzz(k-1)* h(i,j))                             
     &                    /(dzz(k-1)* (h(i,j)+et(i,j)))                   &
! *** NOTE: comment out next line if dens does not include pressure       
     &      +(grav**2)*2.e0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
            gh(i,j,k)=min(gh(i,j,k),.028e0)
          end do
        end do
      end do
! 
      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.e0
          gh(i,j,1)=0.e0
          gh(i,j,kb)=0.e0
        end do
      end do
! 
!    Calculate production of turbulent kinetic energy:
! 
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.25e0*sef                               &
     &                   *((u(i,j,k)-u(i,j,k-1)                           &
     &                      +u(i+1,j,k)-u(i+1,j,k-1))**2                  &
     &                     +(v(i,j,k)-v(i,j,k-1)                          &
     &                      +v(i,j+1,k)-v(i,j+1,k-1))**2)                 &
     &                   /(dzz(k-1)*dh(i,j))**2                           &
!   Add shear due to internal wave field                                  
     &             -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do
! 
!  NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000),
!  depends on ghc the critical number (empirical -6 to -2) to increase mixing.
      ghc=-6.0e0
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.e0
! It is unclear yet if diss. corr. is needed when surf. waves are included.
!           if(gh(i,j,k).lt.0.e0)
!    $        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
!           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)                  &
     &                   /(b1*l(i,j,k)+small)
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))          &
     &                      -(2.e0*dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.e0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)        &
     &                 -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
!     The following section solves the equation:
! 
!       dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
! 
      do j=1,jm
        do i=1,im
          vf(i,j,1)=0.
          vf(i,j,kb)=0.
          ee(i,j,2)=0.e0
          gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
          vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
      end do
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)                                       &
     &                   *(1.e0+e2*((1.e0/abs(z(k)-z(1))                  &
     &                               +1.e0/abs(z(k)-z(kb)))               &
     &                                *l(i,j,k)/(dh(i,j)*kappa))**2)
          end do
        end do
      end do
      do k=3,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))          &
     &                      -(dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)                    &
     &                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
! The following is to counter the problem of the ratio of two
! small numbers (l = q2l/q2) or one number becoming negative.
! Two options are included below. In this application, the second
! option, l  was less noisy when uf or vf is small.
      do k=2,kbm1
        do j=1,jm
          do i=1,im
!           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
!             uf(i,j,k)=small
!             vf(i,j,k)=0.1*dt(i,j)*small
!           endif
          uf(i,j,k)=abs(uf(i,j,k))
          vf(i,j,k)=abs(vf(i,j,k))
          end do
        end do
      end do
! 
! -----------------------------------------------------------------------
! 
!     The following section solves for km and kh:
! 
      coef4=18.e0*a1*a1+9.e0*a1*a2
      coef5=9.e0*a1*a2
! 
!     Note that sm and sh limit to infinity when gh approaches 0.0288:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.e0-6.e0*a1/b1*stf(i,j,k))
            coef2=3.e0*a2*b2/stf(i,j,k)+18.e0*a1*a2
            coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.e0-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.e0-coef5*gh(i,j,k))
          end do
        end do
      end do
! 
!  There are 2 options for kq which, unlike km and kh, was
!  was not derived by Mellor and Yamada but was purely
!  empirical based on neutral boundary layer data.
!  The choice is whether or not it should be subject to
!  the stability factor, sh. Generally, there is not a great
!  difference in output.
      do k=1,kb
        do j=1,jm
          do i=1,im
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:See O2006, Appendix A:                                       !
!                                                                      !
            ZWADMODIFY_C
!           kn(i,j,k)=           l(i,j,k) *sqrt(abs(q2(i,j,k)))
            kn(i,j,k)=(kappa*zsh+l(i,j,k))*sqrt(abs(q2(i,j,k)))
!                                                                      !
! ---------------------------------------------------------------------!
! 
!           kq(i,j,k)=(kn(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
            kq(i,j,k)=(kn(i,j,k)*.20+kq(i,j,k))*.5e0
            km(i,j,k)=(kn(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
            kh(i,j,k)=(kn(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
          end do
        end do
      end do
! cosmetics: make boundr. values as interior
! (even if not used, printout otherwise may show strange values)
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
! 
      return
! 
      end subroutine
! 
! ---------------------------------------------------------------------
! 
      subroutine proft(f,wfsurf,fsurf,nbc,dh)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Solves for vertical diffusion of temperature and    *
! *                salinity using method described by Richmeyer and    *
! *                Morton.                                             *
! *                                                                    *
! *                Irradiance parameters are from Paulson and Simpson. *
! *                                                                    *
! *                See:                                                *
! *                                                                    *
! *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
! *                  Methods for Initial-Value Problems, 2nd edition,  *
! *                  Interscience, New York, 198-201.                  *
! *                                                                    *
! *                Paulson, C. A., and J. Simpson, 1977: Irradiance    *
! *                  measurements in the upper ocean, J. Phys.         *
! *                  Oceanogr., 7, 952-956.                            *
! *                                                                    *
! *                NOTES:                                              *
! *                                                                    *
! *                (1) wfsurf and swrad are negative values when water *
! *                    column is warming or salt is being added.       *
! *                                                                    *
! *                (2) nbc may only be 1 and 3 for salinity.           *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real f(im,jm,kb),wfsurf(im,jm)
      real fsurf(im,jm),dh(im,jm)
      integer nbc
      real rad(im,jm,kb),r(5),ad1(5),ad2(5)
      integer i,j,k,ki
! 
! -----------------------------------------------------------------------
! 
!     Irradiance parameters after Paulson and Simpson:
! 
!       ntp               1      2       3       4       5
!   Jerlov type           i      ia      ib      ii     iii
! 
      data r   /       .58e0,  .62e0,  .67e0,  .77e0,  .78e0 /
      data ad1 /       .35e0,  .60e0,  1.0e0,  1.5e0,  1.4e0 /
      data ad2 /       23.e0,  20.e0,  17.e0,  14.e0,  7.9e0 /
! 
! -----------------------------------------------------------------------
! 
!     Surface boundary condition:
! 
!       nbc   prescribed    prescribed   short wave
!             temperature      flux      penetration
!             or salinity               (temperature
!                                           only)
! 
!        1        no           yes           no
!        2        no           yes           yes
!        3        yes          no            no
!        4        yes          no            yes
! 
!     NOTE that only 1 and 3 are allowed for salinity.
! 
! -----------------------------------------------------------------------
! 
!     The following section solves the equation:
! 
!       dti2*(kh*f')'-f=-fb
! 
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)                             &
     &                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)                               &
     &                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
! 
!     Calculate penetrative radiation. At the bottom any unattenuated
!     radiation is deposited in the bottom layer:
! 
      do k=1,kb
        do j=1,jm
          do i=1,im
            rad(i,j,k)=0.e0
          end do
        end do
      end do
! 
      if(nbc.eq.2.or.nbc.eq.4) then
! 
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=swrad(i,j)                                       &
     &                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp))             &
     &                      +(1.e0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
            end do
          end do
        end do
! 
      endif
! 
      if(nbc.eq.1) then
! 
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do
! 
      else if(nbc.eq.2) then
! 
        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))            &
     &                 /(dz(1)*dh(i,j))                                   &
     &                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do
! 
      else if(nbc.eq.3.or.nbc.eq.4) then
! 
        do j=1,jm
          do i=1,im
            ee(i,j,1)=0.e0
            gg(i,j,1)=fsurf(i,j)
          end do
        end do
! 
      endif
! 
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)                      &
     &                 +dti2*(rad(i,j,k)-rad(i,j,k+1))                    &
     &                   /(dh(i,j)*dz(k)))                                &
     &                 *gg(i,j,k)
          end do
        end do
      end do
! 
!     Bottom adiabatic boundary condition:
! 
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)               &
     &                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))                  &
     &                   /(dh(i,j)*dz(kbm1)))                             &
     &                 /(c(i,j,kbm1)*(1.e0-ee(i,j,kbm2))-1.e0)
        end do
      end do
! 
      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
          f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine profu
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Solves for vertical diffusion of x-momentum using   *
! *                method described by Richmeyer and Morton.           *
! *                                                                    *
! *                See:                                                *
! *                                                                    *
! *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
! *                  Methods for Initial-Value Problems, 2nd edition,  *
! *                  Interscience, New York, 198-201.                  *
! *                                                                    *
! *                NOTE that wusurf has the opposite sign to the wind  *
! *                speed.                                              *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
      real dh(im,jm)
      integer i,j,k,ki
! 
!     The following section solves the equation:
! 
!       dti2*(km*u')'-u=-ub
! 
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do
! 
      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5e0
        end do
      end do
! 
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5e0
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)                              &
     &                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)                                &
     &                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
! 
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))                   &
     &               -uf(i,j,1))                                          &
     &               /(a(i,j,1)-1.e0)
        end do
      end do
! 
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
! 
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i-1,j))                            &
     &              *sqrt(ub(i,j,kbm1)**2                                 &
     &                +(.25e0*(vb(i,j,kbm1)+vb(i,j+1,kbm1)                &
     &                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))            &
     &                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0          &
     &                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do
! 
! 
      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do
! 
      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine profv
! **********************************************************************
!                                                                      *
! * FUNCTION    :  Solves for vertical diffusion of y-momentum using   *
! *                method described by Richmeyer and Morton.           *
! *                                                                    *
! *                See:                                                *
! *                                                                    *
! *                Richtmeyer R.D., and K.W. Morton, 1967. Difference  *
! *                  Methods for Initial-Value Problems, 2nd edition,  *
! *                  Interscience, New York, 198-201.                  *
! *                                                                    *
! *                NOTE that wvsurf has the opposite sign to the wind  *
! *                speed.                                              *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
      real dh(im,jm)
      integer i,j,k,ki
! 
!     The following section solves the equation:
! 
!       dti2*(km*u')'-u=-ub
! 
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do
! 
      do j=2,jm
        do i=2,im
          dh(i,j)=.5e0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do
! 
      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5e0
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)                              &
     &                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)                                &
     &                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do
! 
      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))        &
     &               /(a(i,j,1)-1.e0)
        end do
      end do
! 
      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do
! 
      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i,j-1))                            &
     &              *sqrt((.25e0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)             &
     &                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2   &
     &                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))            &
     &                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0          &
     &                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do
! 
      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do
! 
      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine prxy(label,time,a,im,iskp,jm,jskp,scala)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Writes a horizontal 2-D field.                      *
! *                                                                    *
! *                label ....... label for output                      *
! *                time ........ time (days)                           *
! *                a(im,jm,kb).. array to be printed                   *
! *                iskp ........ skipping interval for i               *
! *                jskp ........ skipping interval for j               *
! *                scala ....... < 0 for floating point numbers output *
! *                              0 for integer output, divisor for a   *
! *                                based on magnitudes of |a| values   *
! *                              > 0 for integer output, divisor for a *
! *                                given by scala                      *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      integer im,jm
      real a(im,jm)
      real time,scala
      integer iskp,jskp
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,jwr,cols
! 
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
! 
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do j=1,jm,jskp
          do i=1,im,iskp
            amx=max(abs(a(i,j)),amx)
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
! 
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
! 
      do ib=1,im,cols*iskp
! 
        ie=ib+(cols-1)*iskp
        if(ie.gt.im) ie=im
! 
        if(scala.ge.0.e0) then
          write(6,3) (i,i=ib,ie,iskp)
    3     format(/,2x,24i5,/)
        else
          write(6,4) (i,i=ib,ie,iskp)
    4     format(/,12i10,/)
        endif
! 
        do j=1,jm,jskp
          jwr=jm+1-j
          if(scala.ge.0.e0) then
            write(6,5) jwr,(nint(a(i,jwr)/scale),i=ib,ie,iskp)
    5       format(1x,i3,24i5)
          else
            write(6,6) jwr,(a(i,jwr),i=ib,ie,iskp)
    6       format(1x,i2,12(e10.2))
          endif
        end do
! 
        write(6,7)
    7   format(//)
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine prxyz(label,time,a,im,iskp,jm,jskp,kb,ko,nko,scala)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Writes horizontal layers of a 3-D field with        *
! *                integers or floating point numbers.                 *
! *                                                                    *
! *                label ....... label for output                      *
! *                time ........ time (days)                           *
! *                a(im,jm,kb).. array to be printed                   *
! *                iskp ........ skipping interval for i               *
! *                jskp ........ skipping interval for j               *
! *                ko .......... 1-D array of k-indices for output     *
! *                nko ......... number of elements in ko              *
! *                scala ....... < 0 for floating point numbers output *
! *                              0 for integer output, divisor for a   *
! *                                based on magnitudes of |a| values   *
! *                              > 0 for integer output, divisor for a *
! *                                given by scala                      *
! *                                                                    *
! *                (NOTE that this combines functions of old prxyz and *
! *                 eprxyz)                                            *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      integer im,jm,kb
      real a(im,jm,kb)
      real time,scala
      integer ko(*)
      integer iskp,jskp,nko
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,jwr,k,iko,cols
! 
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
! 
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do iko=1,nko
          k=ko(iko)
          do j=1,jm,jskp
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
! 
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
! 
      do iko=1,nko
! 
        k=ko(iko)
! 
        write(6,3) k
    3   format(3x,/' Layer k = ',i2)
! 
        do ib=1,im,cols*iskp
! 
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
! 
          if(scala.ge.0.e0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,2x,24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,12i10,/)
          endif
! 
          do j=1,jm,jskp
            jwr=jm+1-j
            if(scala.ge.0.e0) then
              write(6,6) jwr,(nint(a(i,jwr,k)/scale),i=ib,ie,iskp)
    6         format(1x,i3,24i5)
            else
              write(6,7) jwr,(a(i,jwr,k),i=ib,ie,iskp)
    7         format(1x,i2,12(e10.2))
            endif
          end do
! 
          write(6,8)
    8     format(//)
! 
        end do
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine prxz(label,time,a,im,iskp,jm,kb,jo,njo,scala,dt,zz)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
! *                x- or i-direction .                                 *
! *                                                                    *
! *                label ....... label for output                      *
! *                time ........ time (days)                           *
! *                a(im,jm,kb).. array to be printed                   *
! *                iskp ........ skipping interval for i               *
! *                jo .......... 1-D array of j-indices for output     *
! *                njo ......... number of elements in jo              *
! *                scala ....... < 0 for floating point numbers output *
! *                              0 for integer output, divisor for a   *
! *                                based on magnitudes of |a| values   *
! *                              > 0 for integer output, divisor for a *
! *                                given by scala                      *
! *                dt(im,jm) ... total depth                           *
! *                zz(kb) ...... sigma coordinate                      *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      integer im,jm,kb
      real a(im,jm,kb),dt(im,jm),zz(kb)
      real time,scala
      integer jo(*)
      integer iskp,njo
      character label*(*)
      real amx,scale
      integer i,ib,ie,j,k,ijo,cols
! 
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
! 
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do  k=1,kb
          do ijo=1,njo
            j=jo(ijo)
            do i=1,im,iskp
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
! 
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
! 
      do ijo=1,njo
! 
        j=jo(ijo)
! 
        write(6,3) j
    3   format(3x,/' Section j =',i3)
! 
        do ib=1,im,cols*iskp
! 
          ie=ib+(cols-1)*iskp
          if(ie.gt.im) ie=im
! 
          if(scala.ge.0.e0) then
            write(6,4) (i,i=ib,ie,iskp)
    4       format(/,'    i =  ',24i5,/)
          else
            write(6,5) (i,i=ib,ie,iskp)
    5       format(/,'    i =  ',12i10,/)
          endif
! 
          write(6,6) (nint(dt(i,j)),i=ib,ie,iskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
! 
          do k=1,kb
            if(scala.ge.0.e0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),i=ib,ie,iskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),i=ib,ie,iskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
! 
          write(6,9)
    9     format(//)
! 
        end do
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine pryz(label,time,a,im,jm,jskp,kb,io,nio,scala,dt,zz)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
! *                y- or j-direction.                                  *
! *                                                                    *
! *                label ....... label for output                      *
! *                time ........ time (days)                           *
! *                a(im,jm,kb).. array to be printed                   *
! *                jskp ........ skipping interval for j               *
! *                io .......... 1-D array of i-indices for output     *
! *                nio ......... number of elements in io              *
! *                scala ....... < 0 for floating point numbers output *
! *                              0 for integer output, divisor for a   *
! *                                based on magnitudes of |a| values   *
! *                              > 0 for integer output, divisor for a *
! *                                given by scala                      *
! *                dt(im,jm) ... total depth                           *
! *                zz(kb) ...... sigma coordinate                      *
! *                                                                    *
! **********************************************************************
! 
      implicit none
      integer im,jm,kb
      real a(im,jm,kb),dt(im,jm),zz(kb)
      real time,scala
      integer io(*)
      integer jskp,nio
      character label*(*)
      real amx,scale
      integer i,j,jb,je,k,iio,cols
! 
      if(scala.ge.0.e0) then
        cols=24
      else
        cols=12
      endif
! 
      if (scala.lt.0.e0) scale = 1.e0
      if (scala.eq.0.e0) then
        amx=1.e-12
        do  k=1,kb
          do j=1,jm,jskp
            do iio=1,nio
              i=io(iio)
              amx=max(abs(a(i,j,k)),amx)
            end do
          end do
        end do
          scale=10.e0**(int(log10(amx)+100.e0)-103)
        endif
      if(scala.gt.0.e0) scale=scala
! 
      write(6,1) label
    1 format(1x,a40/)
      write(6,2) time,scale
    2 format(' Time = ',f9.4,' days    multiply all values by ',1pe8.2)
! 
      do iio=1,nio
! 
        i=io(iio)
! 
        write(6,3) i
    3   format(3x,/' Section i =',i3)
! 
        do jb=1,jm,cols*jskp
! 
          je=jb+(cols-1)*jskp
          if(je.gt.jm) je=jm
! 
          if(scala.ge.0.e0) then
            write(6,4) (j,j=jb,je,jskp)
    4       format(/,'    j =  ',24i5,/)
          else
            write(6,5) (j,j=jb,je,jskp)
    5       format(/,'    j =  ',12i10,/)
          endif
! 
          write(6,6) (nint(dt(i,j)),j=jb,je,jskp)
    6     format(8x,'d =',24i5.0,/,'     z or zz')
! 
          do k=1,kb
            if(scala.ge.0.e0) then
              write(6,7) k,zz(k),(nint(a(i,j,k)/scale),j=jb,je,jskp)
    7         format(1x,i2,2x,f6.3,24i5)
            else
              write(6,8) k,zz(k),(a(i,j,k),j=jb,je,jskp)
    8         format(1x,i2,2x,f6.3,12(e10.2))
            endif
          end do
! 
          write(6,9)
    9     format(//)
! 
        end do
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine seamount
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Sets up for seamount problem.                       *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real delh,delx,elejmid,elwjmid,ra,vel
      integer i,j,k
! 
!     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
! 
      delh=0.9e0
! 
!     Grid size:
! 
      ZWADMODIFY_C
      delx=8000.e0
! 
!     Radius island or seamount:
! 
      ra=25000.e0
! 
!     Current velocity:
! 
      vel=0.2e0
! 
!     Set up grid dimensions, areas of free surface cells, and
!     Coriolis parameter:
! 
      do j=1,jm
        do i=1,im
! 
!     For constant grid size:
! 
!         dx(i,j)=delx
!         dy(i,j)=delx
! 
!     For variable grid size:
! 
          dx(i,j)=delx-delx*sin(pi*float(i)/float(im))/2.e0
          dy(i,j)=delx-delx*sin(pi*float(j)/float(jm))/2.e0
! 
          cor(i,j)=1.e-4
! 
        end do
      end do
! 
!     Calculate horizontal coordinates of grid points and rotation
!     angle.
! 
!     NOTE that this is introduced solely for the benefit of any post-
!     processing software, and in order to conform with the requirements
!     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
! 
!     There are four horizontal coordinate systems, denoted by the
!     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
!     "e" is an elevation point and "c" is a cell corner), as shown
!     below. In addition, "east_*" is an easting and "north_*" is a
!     northing. Hence the coordinates of the "u" points are given by
!     (east_u,north_u).
! 
!     Also, if the centre point of the cell shown below is at
!     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
!     the coordinates of the western of the two "u" points,
!     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
!     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
!     coordinates of the southwestern corner point of the cell. The
!     southwest corner of the entire grid is at
!     (east_c(1,1),north_c(1,1)).
! 
!                      |              |
!                    --c------v-------c--
!                      |              |
!                      |              |
!                      |              |
!                      |              |
!                      u      e       u
!                      |              |
!                      |              |
!                      |              |
!                      |              |
!                    --c------v-------c--
!                      |              |
! 
! 
!     NOTE that the following calculation of east_c and north_c only
!     works properly for a rectangular grid with east and north aligned
!     with i and j, respectively:
! 
      do j=1,jm
        east_c(1,j)=0.e0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
! 
      do i=1,im
        north_c(i,1)=0.e0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
! 
!     The following works properly for any grid:
! 
!     Elevation points:
! 
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)                          &
     &                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)                       &
     &                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im-1
        east_e(i,jm)                                                      &
     &    =((east_c(i,jm)+east_c(i+1,jm))*3.e0                            &
     &       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
        north_e(i,jm)                                                     &
     &    =((north_c(i,jm)+north_c(i+1,jm))*3.e0                          &
     &       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
      end do
! 
      do j=1,jm-1
        east_e(im,j)                                                      &
     &    =((east_c(im,j)+east_c(im,j+1))*3.e0                            &
     &       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
        north_e(im,j)                                                     &
     &    =((north_c(im,j)+north_c(im,j+1))*3.e0                          &
     &       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
      end do
! 
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)                       &
     &               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)                    &
     &               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
! 
!     u-points:
! 
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
! 
!     v-points:
! 
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
! 
!     rot is the angle (radians, anticlockwise) of the i-axis relative
!     to east, averaged to a cell centre:
! 
!     (NOTE that the following calculation of rot only works properly
!     for this particular rectangular grid)
! 
      do j=1,jm
        do i=1,im
          rot(i,j)=0.e0
        end do
      end do
! 
!     Define depth:
! 
      ZWADMODIFY_C
      hmax=4500.e0   ! tne:!wad: deep open ocean (orig val - not wad-related)
      do i=1,im
        do j=1,jm
! 
      ZWADMODIFY_C
!         h(i,j)=500.*(1.e0-delh                                          &
          h(i,j)=hmax*(1.e0-delh                                          &
     &                          *exp(-((east_c(i,j)                       &
     &                                   -east_c((im+1)/2,j))**2          &
     &                                +(north_c(i,j)                      &
     &                                   -north_c(i,(jm+1)/2))**2)        &
     &                                /ra**2))
          if(h(i,j).lt.1.e0) h(i,j)=1.e0
! 
        end do
      end do
! 
!     Close the north and south boundaries to form a channel:
! 
      do i=1,im
        h(i,1)=1.e0
        h(i,jm)=1.e0
      end do
! 
!     Calculate areas and masks:
! 
      call areas_masks
! 
!     Adjust bottom topography so that cell to cell variations
!     in h do not exceed parameter slmax:
! 
      if(slmax.lt.1.e0) call slpmax
! 
!     Set initial conditions:
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
            ub(i,j,k)=vel*dum(i,j)
          end do
        end do
      end do
! 
!     Initialise uab and vab as necessary
!     (NOTE that these have already been initialised to zero in the
!     main program):
! 
      do j=1,jm
        do i=1,im
          uab(i,j)=vel*dum(i,j)
        end do
      end do
! 
!     Set surface boundary conditions, e_atmos, vflux, wusurf,
!     wvsurf, wtsurf, wssurf and swrad, as necessary
!     (NOTE:
!      1. These have all been initialised to zero in the main program.
!      2. The temperature and salinity of inflowing water must be
!         defined relative to tbias and sbias.):
! 
      do j=1,jm
        do i=1,im
!     No conditions necessary for this problem
! 
! lyo:!wad:!pom2k_bug:tsurf and ssurf were never defined, but should be:
             ZWADMODIFY_C
             tsurf(i,j)=tb(i,j,1)
             ssurf(i,j)=sb(i,j,1)
! 
        end do
      end do
! 
!     Initialise elb, etb, dt and aam2d:
! 
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
! 
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      ZWADMODIFY_C
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
! ---------------------------------------------------------------------!
!                                                                      !
      call dens(sb,tb,rho)
! 
!     Generated horizontally averaged density field (in this
!     application, the initial condition for density is a function
!     of z (the vertical cartesian coordinate) -- when this is not
!     so, make sure that rmean has been area averaged BEFORE transfer
!     to sigma coordinates):
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
! 
!     Set lateral boundary conditions, for use in subroutine bcond
!     (in the seamount problem, the east and west boundaries are open,
!     while the south and north boundaries are closed through the
!     specification of the masks fsm, dum and dvm):
! 
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
! 
      do j=2,jmm1
        uabw(j)=uab(2,j)
        uabe(j)=uab(imm1,j)
! 
!     Set geostrophically conditioned elevations at the boundaries:
! 
        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
      end do
! 
!     Adjust boundary elevations so that they are zero in the middle
!     of the channel:
! 
      elejmid=ele(jmm1/2)
      elwjmid=elw(jmm1/2)
      do j=2,jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j)
        elw(j)=(elw(j)-elwjmid)*fsm(2,j)
      end do
! 
!     Set thermodynamic boundary conditions (for the seamount
!     problem, and other possible applications, lateral thermodynamic
!     boundary conditions are set equal to the initial conditions and
!     are held constant thereafter - users may, of course, create
!     variable boundary conditions):
! 
      do k=1,kbm1
! 
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
! 
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine slpmax
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Limits the maximum of:                              *
! *                                                                    *
! *                  <difference of depths>/<sum of depths>            *
! *                                                                    *
! *                for two adjacent cells. The maximum possible value  *
! *                is unity.                                           *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real mean,del
      integer i,j,loop
! 
      do loop=1,10
! 
!     Sweep right:
! 
        do j=2,jm-1
! 
          do i=2,im-1
            if(fsm(i,j).ne.0.e0.and.fsm(i+1,j).ne.0.e0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.e0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
! 
!    Sweep left:
! 
          do i=im-1,2,-1
            if(fsm(i,j).ne.0.e0.and.fsm(i+1,j).ne.0.e0) then
              if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)).ge.slmax) then
                mean=(h(i+1,j)+h(i,j))/2.e0
                del=sign(slmax,h(i+1,j)-h(i,j))
                h(i+1,j)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
! 
        end do
! 
!   Sweep up:
! 
        do i=2,im-1
! 
          do j=2,jm-1
            if(fsm(i,j).ne.0.e0.and.fsm(i,j+1).ne.0.e0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.e0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
! 
!   Sweep down:
! 
          do j=jm-1,2,-1
            if(fsm(i,j).ne.0.e0.and.fsm(i,j+1).ne.0.e0) then
              if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)).ge.slmax) then
                mean=(h(i,j+1)+h(i,j))/2.e0
                del=sign(slmax,h(i,j+1)-h(i,j))
                h(i,j+1)=mean*(1.e0+del)
                h(i,j)=mean*(1.e0-del)
              endif
            endif
          end do
! 
        end do
! 
      end do
! 
      return
! 
      end subroutine
! 
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff,sw)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates the antidiffusive velocity used to       *
! *                reduce the numerical diffusion associated with the  *
! *                upstream differencing scheme.                       *
! *                                                                    *
! *                This is based on a subroutine of Gianmaria Sannino  *
! *                (Inter-university Computing Consortium, Rome, Italy)*
! *                and Vincenzo Artale (Italian National Agency for    *
! *                New Technology and Environment, Rome, Italy),       *
! *                downloaded from the POM FTP site on 1 Nov. 2001.    *
! *                The calculations have been simplified by removing   *
! *                the shock switch option.                            *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real ff(im,jm,kb)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      real sw
      real mol,abs_1,abs_2
      real value_min,epsilon
      real udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
! 
      parameter (value_min=1.e-9,epsilon=1.0e-14)
! 
!     Apply temperature and salinity mask:
! 
      do k=1,kb
        do i=1,im
          do j=1,jm
            ff(i,j,k)=ff(i,j,k)*fsm(i,j)
          end do
        end do
      end do
! 
!     Recalculate mass fluxes with antidiffusion velocity:
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.                                 &
     &         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.e0
            else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.e0            &
     &              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))                                 &
     &             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.                                 &
     &         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.e0
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.e0             &
     &             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))                                  &
     &            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
! 
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.                                 &
     &         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.e0
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))                                 &
     &             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
      subroutine vertvl(xflux,yflux)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates vertical velocity.                       *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      INCVARDEF
! 
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k
! 
!     Reestablish boundary conditions:
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j))                        &
     &                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do
! 
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1))                        &
     &                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do
! 
!     NOTE that, if one wishes to include freshwater flux, the
!     surface velocity should be set to vflux(i,j). See also
!     change made to 2-D volume conservation equation which
!     calculates elf.
! 
        do j=2,jmm1
          do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
          end do
        end do
! 
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)                                           &
     &                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)                &
     &                        +yflux(i,j+1,k)-yflux(i,j,k))               &
     &                        /(dx(i,j)*dy(i,j))                          &
     &                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do
! 
      return
! 
      end subroutine
! 
! lyo:!wad:                                                             !
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
!     Here are WAD-related subroutines in alphabatical order;          !
!     all names begin with "wad"                                       !
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
      subroutine wadadvt2d(fb,f,ff,nitera,sw,dx,dy,u,v,art,aru,arv,      &
     &    fsm,dum,dvm,aam,dti2)
!                                                                      !
!     2-d version of Smolarkiewicz from /home/lyo/pom/pom2k/           !
!     pom2k.f's subroutine advt2:                                      !
!                                                                      !
!     This version gets rid of common block, i.e. the routine is now   !
!     self-contained, except that im & jm need to be specified below   !
!     if "include 'grid'" is NOT used; also, 2-d calculation, i.e,     !
!     solve for  D:                                                    !
!                                                                      !
!     d(D)/dt + d(UD)/dx + d(VD)/dy = Diffusion_of_(D)                 !
!                                                                      !
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Integrates conservative scalar equations.           *
! *                                                                    *
! *                This is a first-order upstream scheme, which        *
! *                reduces implicit diffusion using the Smolarkiewicz  *
! *                iterative upstream scheme with an antidiffusive     *
! *                velocity.                                           *
! *                                                                    *
! *                It is based on the subroutines of Gianmaria Sannino *
! *                (Inter-university Computing Consortium, Rome, Italy)*
! *                and Vincenzo Artale (Italian National Agency for    *
! *                New Technology and Environment, Rome, Italy),       *
! *                downloaded from the POM FTP site on 1 Nov. 2001.    *
! *                The calculations have been simplified by removing   *
! *                the shock switch option. It should be noted that    *
! *                this implementation does not include cross-terms    *
! *                which are in the original formulation.              *
! *                                                                    *
! *                fb,f,fclim,ff . as used in subroutine advt1         *
! *                xflux,yflux ... working arrays used to save memory  *
! *                nitera ........ number of iterations. This should   *
! *                                be in the range 1 - 4. 1 is         *
! *                                standard upstream differencing;     *
! *                                3 adds 50% CPU time to POM.         *
! *                sw ............ smoothing parameter. This should    *
! *                                preferably be 1, but 0 < sw < 1     *
! *                                gives smoother solutions with less  *
! *                                overshoot when nitera > 1.          *
! *                                                                    *
! *                Reference:                                          *
! *                                                                    *
! *                Smolarkiewicz, P.K.; A fully multidimensional       *
! *                  positive definite advection transport algorithm   *
! *                  with small implicit diffusion, Journal of         *
! *                  Computational Physics, 54, 325-362, 1984.         *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      real sw,dti2
      integer nitera,im,jm,kb  ! lyo:!wad:kb not required but defined
                               !         in "grid" below
! 
!     PARAMETER (IM=65,JM=49)
!     PARAMETER (IM=131,JM=99)
      include 'grid'
! 
      real fb(im,jm),f(im,jm),ff(im,jm)
      real dx(im,jm),dy(im,jm),u(im,jm),v(im,jm),aam(im,jm)
      real fsm(im,jm),dum(im,jm),dvm(im,jm)
      real art(im,jm),aru(im,jm),arv(im,jm)

!     integer ix,jx
!     PARAMETER (IX=62,JX=18)
      real xflux(im,jm),yflux(im,jm)
      real fbmem(im,jm)
      real xmassflux(im,jm),ymassflux(im,jm)
      integer i,j,k,itera,imm1,jmm1
! 
!     Calculate horizontal mass fluxes:
! 
      imm1=im-1; jmm1=jm-1
! 
        do j=1,jm
          do i=1,im
            xmassflux(i,j)=0.e0
            ymassflux(i,j)=0.e0
          end do
        end do
! 
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j)=0.5e0*(dy(i-1,j)+dy(i,j))*u(i,j)
          end do
        end do
! 
        do j=2,jm
          do i=2,imm1
            ymassflux(i,j)=0.5e0*(dx(i,j-1)+dx(i,j))*v(i,j)
          end do
        end do
! 
        do j=1,jm
          do i=1,im
            fbmem(i,j)=fb(i,j)
          end do
        end do
! 
!     Start Smolarkiewicz scheme:
! 
      do itera=1,nitera
! 
!     Upwind advection scheme:
! 
          do j=2,jm
            do i=2,im
              xflux(i,j)=0.5e0                                            &
     &                      *((xmassflux(i,j)+abs(xmassflux(i,j)))        &
     &                        *fbmem(i-1,j)+                              &
     &                        (xmassflux(i,j)-abs(xmassflux(i,j)))        &
     &                        *fbmem(i,j))
! 
              yflux(i,j)=0.5e0                                            &
     &                      *((ymassflux(i,j)+abs(ymassflux(i,j)))        &
     &                        *fbmem(i,j-1)+                              &
     &                        (ymassflux(i,j)-abs(ymassflux(i,j)))        &
     &                        *fbmem(i,j))
            end do
          end do
! 
!     Add net advective fluxes and step forward in time:
! 
          do j=2,jmm1
            do i=2,imm1
              ff(i,j)=xflux(i+1,j)-xflux(i,j)                             &
     &                 +yflux(i,j+1)-yflux(i,j)
              ff(i,j)=(fbmem(i,j)*art(i,j)                                &
     &                   -dti2*ff(i,j))/(art(i,j))
            end do
          end do
! 
!     Calculate antidiffusion velocity:
! 
      call wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2)
!     call wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2,im,jm)
! 
        do j=1,jm
          do i=1,im
              fbmem(i,j)=ff(i,j)
          end do

        end do
! 
!     End of Smolarkiewicz scheme
! 
      end do
! 
!     Add horizontal diffusive fluxes:
! 
        do j=2,jm
          do i=2,im
            xmassflux(i,j)=0.5e0*(aam(i,j)+aam(i-1,j))
            ymassflux(i,j)=0.5e0*(aam(i,j)+aam(i,j-1))
          end do
        end do
! 
        do j=2,jm
          do i=2,im
           xflux(i,j)=-xmassflux(i,j)                                     &
     &                   *(fb(i,j)-fb(i-1,j))*dum(i,j)                    &
     &                   *(dy(i,j)+dy(i-1,j))/(dx(i,j)+dx(i-1,j))
           yflux(i,j)=-ymassflux(i,j)                                     &
     &                   *(fb(i,j)-fb(i,j-1))*dvm(i,j)                    &
     &                   *(dx(i,j)+dx(i,j-1))/(dy(i,j)+dy(i,j-1))
          end do
        end do
! 
!     Add net horizontal fluxes and step forward in time:
! 
        do j=2,jmm1
          do i=2,imm1
            ff(i,j)=ff(i,j)-dti2*(xflux(i+1,j)-xflux(i,j)                 &
     &                               +yflux(i,j+1)-yflux(i,j))            &
     &                           /(art(i,j))
          end do
        end do
! 
      return
! 
      end subroutine
! 
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
      subroutine wadout
!                                                                      !
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  WAD Outputs (temporary subroutine to be merged with *
! *                the rest of the code later)                         *
! *                                                                    *
! *                Edit as approriate.                                 *
! *                                                                    *
! **********************************************************************
!                                                                      !
      implicit none
! 
      integer i,j,k
! 
      INCVARDEF
! 
!     2-D horizontal fields:
! 
          call prxy('wetmask; =0 are land or dry, =1 are wet ',           &
     &              time,wetmask,im,iskp,jm,jskp,1.e0)
! 
          do j=1,jm; do i=1,im
             tps(i,j)=fsm(i,j)-wetmask(i,j)
             enddo; enddo
          call prxy('fsm-wetmask; WAD: =1 are dry cells      ',           &
     &              time,tps,im,iskp,jm,jskp,1.e0)
! 
          do j=1,jm; do i=1,im
             tps(i,j)=(elb(i,j)+hhi)*wetmask(i,j)
             enddo; enddo
          call prxy('elb+hhi (m; w.r.t MSL)                  ',           &
     &              time,tps,im,iskp,jm,jskp,0.e0)
! 
          do j=1,jm; do i=1,im
             tps(i,j)=d(i,j)*fsm(i,j)
             enddo; enddo
          call prxy('Total Water Depth (m)                   ',           &
     &              time,tps,im,iskp,jm,jskp,0.e0)
! 
! tne: -------------- save for matlab plot
! 
       write(88,'(2i5,4f8.3)')im-2,jm-2,time*24.,99.9,99.9,99.9
      do j=2,jm-1; do i=2,im-1
       write(88,'(2i5,4f8.3)')i,j,tps(i,j),uab(i,j),vab(i,j),d(i,j)
      enddo; enddo
! 
      return
      end subroutine
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
      subroutine wadh
!                                                                      !
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Sets up WAD definitions of topography [see Fig.1 of *
! *                O2006]                                              *
! *                                                                    *
! **********************************************************************
!                                                                      !
!     Input (through common block) is "h" defined as follows:          !
!                                                                      !
!     "h" is wrt MSL (Mean Sea Level) at which level h=0, and h<0 for  !
!     depth below the MSL, and h>0 for "land" above the MSL (i.e. for  !
!     marshes, wetlands, hills & mountains etc.)                       !
!                                                                      !
!     Also, nwad, hhi, slmax                                           !
!                                                                      !
!     Outputs:                                                         !
!                                                                      !
!     h --> (i.e. becomes) h_old_pom + hhi (which can be zero)         !
!       i.e. "h" is wrt to some high land-datum (Absolute Land Boundary!
!       or ALB), "hhi" meters above the MSL. All cells other than the  !
!       ALB cells are either always wet or can potentially become wet. !
!                                                                      !
!     Also, fsm, dum, dvm, wetmask, cell areas etc, and hmax           !
!       see more detailed definitions inside the subroutine            !
!                                                                      !
! ---------------------------------------------------------------------!
      implicit none
!                                                                      !
      INCVARDEF
!                                                                      !
      integer npos,nneg
      real hwatmin
      integer i,j
!                                                                      !
! ---------------------------------------------------------------------!
!     Do a rudimentary check on "h:"                                   !
!                                                                      !
      npos=0; nneg=0
      do j=1,jm; do i=1,im
      if (h(i,j).ge.0.0) then
      npos=+1
      else
      nneg=-1
      endif
      if (npos*nneg .lt. 0) go to 101
      enddo; enddo
!                                                                      !
      write(6,'('' Stopped in subr. wadh; incorrect defn of h:'')')
      write(6,'('' h is one-sign only; see comments in subr.wadh'',/)')
      write(6,'('' npos,nneg ='',2i5,/)') npos,nneg
      call prxy('Undisturbed water depth, h',time,h,im,1,jm,1,0.0)
      stop
!                                                                      !
 101  continue
!                                                                      !
      hkeep(:,:)=h(:,:) ! Keep original for plots etc
!                                                                      !
!     HMIN & HC (all +ve) are defined as:                              !
!                                                                      !
!     h = hmin ---> absolute land area [never flooded], FSM=0.0;       !
!     h >= hc  ---> water depth wrt High water level, FSM=1.0;         !
!                                                                      !
!     Note that hc is defined in MAIN --                               !
!     hmin=0.01; hc=0.05; hco=hc       !hco is used to avoid round-
      hmin=0.01; hco=hc                ! hco is used to avoid round-
      hc=hc*1.0001                     ! off in initial elevation
!     hc=hc*1.01                       !tne:!wad:- using a larger
!                 multiplying factor, "1.01" instead of "1.0001",      !
!                 can help eliminate singular wet spots due to roundoff!
!                                                                      !
      write(6,'('' Outputs from subr. wadh ---------------------'',/)')
      write(6,'('' hmin,hc,hco ='',1p3e13.4,/)') hmin,hc,hco
!                                                                      !
!     Reverse sign of "h", and shift reference if hhi.ne.0, see below..!
!     (Note that hhi is defined in MAIN)                               !
      h(:,:)=-h(:,:)+hhi
!                                                                      !
      if (nwad.eq.1) then
!                                                                      !
!     If nwad=1, then hhi.ne.0.0, and line "h(:,:)=-h(:,:)+hhi" above  !
!     already shifts the reference level from MSL to ALB, the Absolute !
!     Land Boundary, according to O2006, Fig.1:                        !
!                                                                      !
      do j=1,jm; do i=1,im
      FSM(I,J)=1.; DUM(I,J)=1.; DVM(I,J)=1.
      IF(h(I,J).LT.0.0)THEN ! Define absolute land areas (never flood):
                            ! temporarily set h=-1.0 (some -ve number)
      h(I,J)=-1.0; FSM(I,J)=0.; DUM(I,J)=0.; DVM(I,J)=0.
      ENDIF
      enddo; enddo
      DO J=1,JM-1; DO I=1,IM
      IF(FSM(I,J).EQ.0..AND.FSM(I,J+1).NE.0.)DVM(I,J+1)=0.
      enddo; enddo
      DO J=1,JM;   DO I=1,IM-1
      IF(FSM(I,J).EQ.0..AND.FSM(I+1,J).NE.0.)DUM(I+1,J)=0.
      enddo; enddo
!                                                                      !
!     The above yields the following definitions for h:                !
!                                                                      !
!     h=-1, absolute land areas [never flooded], FSM=0.0, always       !
!     h>=0, wet but potentially dry areas,  FSM also=1.0, always       !
!                                                                      !
!     Note that slpmax touches fsm=1 points only - wet or potentially  !
!     WAD cells; it changes "h" so a wet cell might become dry; so it  !
!     must be called here before wetmask is defined:                   !
!                                                                      !
!     Adjust bottom topography so that cell to cell variations         !
!     in h do not exceed parameter slmax:                              !
!                                                                      !
      if(slmax.lt.1.e0) call slpmax
!                                                                      !
!     We now define a wetmask, assuming that at t=0, the free surface  !
!     is at the MSL (i.e. elevation=0 wrt MSL).  If this run does not  !
!     begin with time=0 (or if other special free-surface distribution !
!     is desired), then wetmask needs to be redefined or read in e.g.  !
!     from a RESTART file, as done later in MAIN if NREAD.ne.0         !
!                                                                      !
      wetmask(:,:)=fsm(:,:)
      do j=1,jm; do i=1,im
      if (h(i,j).lt.hhi) wetmask(i,j)=0.0
      enddo; enddo
!                                                                      !
!     Finalize definitions of "h":                                     !
!                                                                      !
!     h=hmin, absolute land areas [never flooded], FSM=0.0, always     !
!     h>=hc,  wet but potentially dry areas,  FSM also=1.0, always     !
!                                   ...  but wetmask can be 0 or 1     !
!                                                                      !
      do j=1,jm; do i=1,im
      if (h(i,j).lt.0.0) then
      h(i,j)=hmin
      else
      h(i,j)=h(i,j)+hc
      endif
      enddo; enddo
!                                                                      !
      call areas_masks ! Calc. cells' areas etc but do not alter fsm...
                       ! if nwad=1
!                                                                      !
      else    ! nwad.eq.0
!                                                                      !
      DO J=1,JM; DO I=1,IM
      IF (h(I,J).LE.0.0) h(I,J)=hmin
      ENDDO; ENDDO
!                                                                      !
!     Set minimum water depth to hwatmin:                              !
!     hwatmin=min water depth, must be >1 to match subr.areas_masks'   !
!     definition of land, and large enough to avoid grid cells from    !
!     becoming dry                                                     !
!                                                                      !
      hwatmin=10.0
      DO J=1,JM; DO I=1,IM
      IF (h(I,J).GT.hmin.and.h(I,J).LT.hwatmin) h(I,J)=hwatmin
      ENDDO; ENDDO
!                                                                      !
      call areas_masks      ! Calc. cells' areas etc & fsm if nwad=0
!                                                                      !
!     Adjust bottom topography so that cell to cell variations         !
!     in h do not exceed parameter slmax:                              !
!                                                                      !
      if(slmax.lt.1.e0) call slpmax
!                                                                      !
      wetmask(:,:)=fsm(:,:) ! wetmask is same as fsm if nwad=0
!                                                                      !
      endif   ! if (nwad.eq.1) then...else...
!                                                                      !
!     Print to check:                                                  !
!                                                                      !
      CALL PRXY(' Input h ',TIME, hkeep(1,1),IM,1,JM,1,0.0)
      CALL PRXY(' h after wadh ',TIME, h(1,1),IM,1,JM,1,0.0)
      CALL PRXY(' FSM ',TIME, FSM(1,1),IM,1,JM,1,1.0)
      CALL PRXY(' WETMASK ',TIME, WETMASK(1,1),IM,1,JM,1,1.0)
      tps(:,:)=FSM(:,:)-WETMASK(:,:)
      CALL PRXY('FSM-WETMASK',TIME,tps(1,1),IM,1,JM,1,1.0)
!                                                                      !
!     Double-check masks for nwad=0:                                   !
!                                                                      !
      if (nwad.eq.0) then
      do j=1,jm; do i=1,im
      if (fsm(i,j).ne.wetmask(i,j)) then
      write(6,'('' Stopped, nwad ='',i4)') nwad
      write(6,'('' fsm.ne.wetmask @ (i,j) ='',2i4,/)') i,j
      stop
      endif
      enddo; enddo
      endif
!                                                                      !
      return
      end subroutine
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
      subroutine wadseamount()
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Sets up for seamount problem w/wad.                 *
! *                                                                    *
! **********************************************************************
! 
      implicit none 
! 
      INCVARDEF 
! 
      real delh,delx,elejmid,elwjmid,ra,vel,hland  ! lyo:!wad:define hland
      integer nvar                                 ! lyo:!wad:define nvar
      integer i,j,k
! 
!     Set delh > 1.0 for an island or delh < 1.0 for a seamount:
! 
      delh=1.15e0   ! tne:!wad: island (for wad) =0.9e0 for orig.seamount
! 
!     Grid size:
! 
      delx=8000.e0
! 
! lyo:!wad:Define variable or uniform grid option:
      nvar=1  ! =1 for variable grid; =0 otherwise
! lyo:!wad:Special test case for smaller delx=4km:
      if ((nvar.eq.0).or.(im.eq.131.and.jm.eq.99)) delx=delx*0.5
! 
!     Radius island or seamount:
! 
      ra=50000.e0   ! tne:!wad:
! 
!     Current velocity:
! 
      vel=0.2e0
! 
!     Set up grid dimensions, areas of free surface cells, and
!     Coriolis parameter:
! 
      do j=1,jm
        do i=1,im
! 
!     For constant grid size:
! 
!         dx(i,j)=delx
!         dy(i,j)=delx
! 
!     For variable grid size:
! 
! lyo:!wad:variable or uniform grid--
          dx(i,j)=delx-float(nvar)*delx*sin(pi*float(i)/float(im))/2.e0
          dy(i,j)=delx-float(nvar)*delx*sin(pi*float(j)/float(jm))/2.e0
! 
          cor(i,j)=1.e-4
! 
        end do
      end do
! 
!     Calculate horizontal coordinates of grid points and rotation
!     angle.
! 
!     NOTE that this is introduced solely for the benefit of any post-
!     processing software, and in order to conform with the requirements
!     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
! 
!     There are four horizontal coordinate systems, denoted by the
!     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
!     "e" is an elevation point and "c" is a cell corner), as shown
!     below. In addition, "east_*" is an easting and "north_*" is a
!     northing. Hence the coordinates of the "u" points are given by
!     (east_u,north_u).
! 
!     Also, if the centre point of the cell shown below is at
!     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
!     the coordinates of the western of the two "u" points,
!     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
!     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
!     coordinates of the southwestern corner point of the cell. The
!     southwest corner of the entire grid is at
!     (east_c(1,1),north_c(1,1)).
! 
!                      |              |
!                    --c------v-------c--
!                      |              |
!                      |              |
!                      |              |
!                      |              |
!                      u      e       u
!                      |              |
!                      |              |
!                      |              |
!                      |              |
!                    --c------v-------c--
!                      |              |
! 
! 
!     NOTE that the following calculation of east_c and north_c only
!     works properly for a rectangular grid with east and north aligned
!     with i and j, respectively:
! 
      do j=1,jm
        east_c(1,j)=0.e0
        do i=2,im
          east_c(i,j)=east_c(i-1,j)+dx(i-1,j)
        end do
      end do
! 
      do i=1,im
        north_c(i,1)=0.e0
        do j=2,jm
          north_c(i,j)=north_c(i,j-1)+dy(i,j-1)
        end do
      end do
! 
!     The following works properly for any grid:
! 
!     Elevation points:
! 
      do j=1,jm-1
        do i=1,im-1
          east_e(i,j)=(east_c(i,j)+east_c(i+1,j)                          &
     &                  +east_c(i,j+1)+east_c(i+1,j+1))/4.e0
          north_e(i,j)=(north_c(i,j)+north_c(i+1,j)                       &
     &                   +north_c(i,j+1)+north_c(i+1,j+1))/4.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im-1
        east_e(i,jm)                                                      &
     &    =((east_c(i,jm)+east_c(i+1,jm))*3.e0                            &
     &       -east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0
        north_e(i,jm)                                                     &
     &    =((north_c(i,jm)+north_c(i+1,jm))*3.e0                          &
     &       -north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0
      end do
! 
      do j=1,jm-1
        east_e(im,j)                                                      &
     &    =((east_c(im,j)+east_c(im,j+1))*3.e0                            &
     &       -east_c(im-1,j)-east_c(im-1,j+1))/4.e0
        north_e(im,j)                                                     &
     &    =((north_c(im,j)+north_c(im,j+1))*3.e0                          &
     &       -north_c(im-1,j)-north_c(im-1,j+1))/4.e0
      end do
! 
      east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)                       &
     &               -(east_e(im-2,jm)+east_e(im,jm-2))/2.e0
      north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)                    &
     &               -(north_e(im-2,jm)+north_e(im,jm-2))/2.e0
! 
!     u-points:
! 
      do j=1,jm-1
        do i=1,im
          east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0
          north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do i=1,im
        east_u(i,jm)=(east_c(i,jm)*3.e0-east_c(i,jm-1))/2.e0
        north_u(i,jm)=(north_c(i,jm)*3.e0-north_c(i,jm-1))/2.e0
      end do
! 
!     v-points:
! 
      do j=1,jm
        do i=1,im-1
          east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0
          north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0
        end do
      end do
! 
!     Extrapolate ends:
! 
      do j=1,jm
        east_v(im,j)=(east_c(im,j)*3.e0-east_c(im-1,j))/2.e0
        north_v(im,j)=(north_c(im,j)*3.e0-north_c(im-1,j))/2.e0
      end do
! 
!     rot is the angle (radians, anticlockwise) of the i-axis relative
!     to east, averaged to a cell centre:
! 
!     (NOTE that the following calculation of rot only works properly
!     for this particular rectangular grid)
! 
      do j=1,jm
        do i=1,im
          rot(i,j)=0.e0
        end do
      end do
! 
!     Define depth:
! 
! tne:!wad:!lyo:!wad:
!     Note that unlike original seamount, h<0 is now water below MSL, and
!     h>0 is above the MSL and is either land (if h>hland, see below)
!     or is potential wet-and-dry region
! 
      hmax=50.e0; hland=5.e0 ! note all pnts are potential WAD if we set
                             ! hland >= -hmax*(1.e0-delh) (=7.5m);
                             ! i.e. "if" below is NOT satisfied
! 
      do i=1,im
        do j=1,jm
! 
          h(i,j)=-hmax*(1.e0-delh                                         &
     &                          *exp(-((east_c(i,j)                       &
     &                                   -east_c((im+1)/2,j))**2          &
     &                                +(north_c(i,j)                      &
     &                                   -north_c(i,(jm+1)/2))**2)        &
     &                                /ra**2))
! lyo:!wad:Make region near top of seamount absolute land region:
          if(h(i,j).gt.hland) h(i,j)=hhi+1.e0
! 
        end do
      end do
! 
!     Close the north and south boundaries to form a channel:
! 
      do i=1,im
        h(i,1)=hhi+1.e0  ! tne:!wad:
        h(i,jm)=hhi+1.e0 ! tne:!wad:
      end do
! 
!     Calculate areas and masks:
! 
! lyo:!wad:Define "h" consistent with WAD and calculate wetmask, cell areas
!         and fsm etc.
!     call areas_masks
      call wadh
! 
!     Adjust bottom topography so that cell to cell variations
!     in h do not exceed parameter slmax:
! 
!     if(slmax.lt.1.e0) call slpmax  !lyo:wad:now done in wadh above
! 
!     Set initial conditions:
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
! lyo:wad:
!           tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
             if (h(i,j).gt.hhi) then ! lyo:wad:stratified water below MSL:
                tb(i,j,k)=5.e0+15.e0*exp(zz(k)*(h(i,j)-hhi)/1000.e0)
             else                    ! lyo:wad:well-mixed water above MSL:
                tb(i,j,k)=5.e0+15.e0
             endif
            tb(i,j,k)=tb(i,j,k)-tbias
! 
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
            ub(i,j,k)=vel*dum(i,j)
          end do
        end do
      end do
! 
!     Initialise uab and vab as necessary
!     (NOTE that these have already been initialised to zero in the
!     main program):
! 
      do j=1,jm
        do i=1,im
          uab(i,j)=vel*dum(i,j)
        end do
      end do
! 
!     Set surface boundary conditions, e_atmos, vflux, wusurf,
!     wvsurf, wtsurf, wssurf and swrad, as necessary
!     (NOTE:
!      1. These have all been initialised to zero in the main program.
!      2. The temperature and salinity of inflowing water must be
!         defined relative to tbias and sbias.):
! 
      do j=1,jm
        do i=1,im
!     No conditions necessary for this problem
! 
! lyo:!wad:!pom2k_bug:tsurf and ssurf were never defined, but should be:
             tsurf(i,j)=tb(i,j,1)
             ssurf(i,j)=sb(i,j,1)
! 
        end do
      end do
! 
!     Initialise elb, etb, dt and aam2d:
! 
      do j=1,jm
        do i=1,im
          elb(i,j)=(-e_atmos(i,j)-hhi)*fsm(i,j)      ! lyo:!wad:
! lyo:!wad:note:The following "if" is not satisfied if nwad=0 since     !
!     then fsm=wetmask at all cells                                    !
!     if (fsm(i,j).ne.0.0.and.wetmask(i,j).eq.0.0) ! dry but not land
      if (fsm(i,j).ne.wetmask(i,j))                ! dry but not land     &
     &elb(i,j)=-h(i,j)+hco ! note slightly smaller hco<hc to ensure that
                           ! initially, cell is wet with elb+h=hco < hc
          etb(i,j)=elb(i,j)               ! lyo:!wad:
          dt(i,j)=h(i,j)+elb(i,j)         ! lyo:!wad:
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad:Make sure that initial wetmask remains compatible with       !
!     initial elb.  Also define marsh evaporation, +ve upward or       !
!     evaporate.  Also initialize wriv for rivers.                     !
! lyo:!wad:note:wriv is added here but the implementation following     !
!     Oey [1996; JPO, 26, 145-175; but see p.154 in particular] is not !
!     complete in this code.  E-mail me if want to implement it.       !
!                                                                      !
      if (nwad.eq.1) then
      do j=1,jm; do i=1,im
      wetmask(i,j)=fsm(i,j)
      if ((h(i,j)+elb(i,j)).le.hc) wetmask(i,j)=0.0
      enddo; enddo
      endif
!                                                                      !
      do j=1,jm; do i=1,im
      wmarsh(i,j)=0.0; wriv(i,j)=0.0
      if (fsm(i,j).ne.wetmask(i,j)) then           ! dry but not land
!     2 examples of wmarsh:
! --   wmarsh(i,j)=fsm(i,j)*10./art(i,j)
! --   wmarsh(i,j)=fsm(i,j)*(4.e-7)  !Gill's [1982] value
      endif
      enddo; enddo
!                                                                      !
! ---------------------------------------------------------------------!
! lyo:!wad: Set up pdens before 1st call dens; used also in profq:      !
      do k=1,kbm1; do j=1,jm; do i=1,im
         pdens(i,j,k)=grav*rhoref*(-zz(k)*max(h(i,j)-hhi,0.e0))*1.e-5
      enddo; enddo; enddo
! ---------------------------------------------------------------------!
!                                                                      !
      call dens(sb,tb,rho)
! 
!     Generated horizontally averaged density field (in this
!     application, the initial condition for density is a function
!     of z (the vertical cartesian coordinate) -- when this is not
!     so, make sure that rmean has been area averaged BEFORE transfer
!     to sigma coordinates):
! 
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
! 
!     Set lateral boundary conditions, for use in subroutine bcond
!     (in the seamount problem, the east and west boundaries are open,
!     while the south and north boundaries are closed through the
!     specification of the masks fsm, dum and dvm):
! 
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
! 
      do j=2,jmm1
        uabw(j)=uab(2,j)
        uabe(j)=uab(imm1,j)
! 
!     Set geostrophically conditioned elevations at the boundaries:
! 
! lyo:!wad:Note keep (temporarily) all el* defined wrt MSL - ele, elw,
!         eln & els were initialized to be =0 in MAIN; then adjust
!         them later (below) with hhi:
! 
        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
      end do
! 
!     Adjust boundary elevations so that they are zero in the middle
!     of the channel:
! 
      elejmid=ele(jmm1/2)
      elwjmid=elw(jmm1/2)
      do j=2,jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j)
        elw(j)=(elw(j)-elwjmid)*fsm(2,j)
      end do
! 
! tne:!wad:!lyo:wad:
      do i=1,im
         eln(i)=(0.-hhi)*fsm(i,jm)
         els(i)=(0.-hhi)*fsm(i,2)
      enddo
      do j=1,jm
         ele(j)=(ele(j)-hhi)*fsm(im,j)
         elw(j)=(elw(j)-hhi)*fsm(2,j)
      enddo
! 
!     Set thermodynamic boundary conditions (for the seamount
!     problem, and other possible applications, lateral thermodynamic
!     boundary conditions are set equal to the initial conditions and
!     are held constant thereafter - users may, of course, create
!     variable boundary conditions):
! 
      do k=1,kbm1
! 
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
! 
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
! 
      end do
! 
      return
! 
      end subroutine
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
      subroutine wadsmoladif(xmassflux,ymassflux,ff,sw,fsm,aru,arv,dti2) 
!                                                                      !
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Calculates the antidiffusive velocity used to       *
! *                reduce the numerical diffusion associated with the  *
! *                upstream differencing scheme.                       *
! *                                                                    *
! *                This is based on a subroutine of Gianmaria Sannino  *
! *                (Inter-university Computing Consortium, Rome, Italy)*
! *                and Vincenzo Artale (Italian National Agency for    *
! *                New Technology and Environment, Rome, Italy),       *
! *                downloaded from the POM FTP site on 1 Nov. 2001.    *
! *                The calculations have been simplified by removing   *
! *                the shock switch option.                            *
! *                                                                    *
! **********************************************************************
! 
      implicit none
! 
      integer im,jm,kb  ! lyo:!wad:kb not required but defined
                        !         in "grid" below
! 
!     PARAMETER (IM=65,JM=49)
!     PARAMETER (IM=131,JM=99)
      include 'grid'
! 
      real ff(im,jm)
      real xmassflux(im,jm),ymassflux(im,jm)
      real fsm(im,jm),aru(im,jm),arv(im,jm)
      real sw,dti2
      real mol,abs_1,abs_2
      real value_min,epsilon
      real udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k,imm1,jmm1
! 
      parameter (value_min=1.e-9,epsilon=1.0e-14)
! 
! 
      imm1=im-1; jmm1=jm-1
! 
!     Apply temperature and salinity mask:
! 
        do i=1,im
          do j=1,jm
            ff(i,j)=ff(i,j)*fsm(i,j)
          end do
        end do
! 
!     Recalculate mass fluxes with antidiffusion velocity:
! 
        do j=2,jmm1
          do i=2,im
            if(ff(i,j).lt.value_min.or.                                   &
     &         ff(i-1,j).lt.value_min) then
              xmassflux(i,j)=0.e0
            else
              udx=abs(xmassflux(i,j))
              u2dt=dti2*xmassflux(i,j)*xmassflux(i,j)                     &
     &              /(aru(i,j))
              mol=(ff(i,j)-ff(i-1,j))                                     &
     &             /(ff(i-1,j)+ff(i,j)+epsilon)
              xmassflux(i,j)=(udx-u2dt)*mol*sw

              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j)=0.e0
            end if
          end do
        end do
! 
        do j=2,jm
          do i=2,imm1
            if(ff(i,j).lt.value_min.or.                                   &
     &         ff(i,j-1).lt.value_min) then
              ymassflux(i,j)=0.e0
            else
             vdy=abs(ymassflux(i,j))
             v2dt=dti2*ymassflux(i,j)*ymassflux(i,j)                      &
     &             /(arv(i,j))
             mol=(ff(i,j)-ff(i,j-1))                                      &
     &            /(ff(i,j-1)+ff(i,j)+epsilon)
             ymassflux(i,j)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j)=0.e0
            end if
          end do
        end do
! 
      return
! 
      end subroutine
!                                                                      !
! lyo:!wad:END of WAD-related subroutines.                              !
!                                                                      !
! ---------------------------------------------------------------------!
!                                                                      !
!      include 'pom2008.n'                                     ! *netCDF*
! 
!     End of source code
! 
! -----------------------------------------------------------------------
