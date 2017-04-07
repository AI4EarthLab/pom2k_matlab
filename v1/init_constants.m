im=7; jm=5; kb=6;
% im=21; jm=21; kb=6;
% im=11; jm=9; kb=6;
% im=65; jm=49; kb=21;
imm1=im-1; imm2=im-2; jmm1=jm-1; jmm2=jm-2; kbm1=kb-1; kbm2=kb-2;

alpha          =0.0;dte            =0.0;dti            =0.0;dti2           =0.0;
grav           =0.0;hmax           =0.0;kappa          =0.0;pi             =0.0;
ramp           =0.0;rfe            =0.0;rfn            =0.0;rfs            =0.0;
rfw            =0.0;rhoref         =0.0;sbias          =0.0;slmax          =0.0;
small          =0.0;tbias          =0.0;time           =0.0;tprni          =0.0;
umol           =0.0;vmaxl 		    =0.0;

iint           =0;iprint         =0;iskp           =0;jskp         =0;
kl1            =0;kl2            =0;mode           =0;ntp			=0;

%     Integers defining the number of logarithmic layers at the
%     surface and bottom (used by subroutine depth). The number of
%     logarithmic layers are kl1-2 at the surface and kb-kl2-1
%     at the bottom. For no log portions, set kl1=2 and kl2=kb-1:
kl1=6;
kl2=kb-2;

source='pom2k  2006-05-03';

small=1.e-9;          							% Small value

pi=atan(1.0)*4.0;    							% PI

title='Run 1'; 									% run's title

netcdf_file='pom2k.nc';  						% netCDF output file

%     Problem number:
%     iproblem      problem      initialisation
%                    type          subroutine
%         1        seamount       seamount
%         2        conservation   box
%                  box
%         3        IC from file   file2ic
iproblem=1;

%       mode                     description
%        2        2-D calculation (bottom stress calculated in advave)
%        3        3-D calculation (bottom stress calculated in profu,v)
%        4        3-D calculation with t and s held fixed
mode=3;

%     Advection scheme:
%      nadv     Advection scheme
%        1       Centred scheme, as originally provide in POM
%        2       Smolarkiewicz iterative upstream scheme, based on
%                subroutines provided by Gianmaria Sannino and Vincenzo
%                Artale
nadv=1;

%     Constants for Smolarkiewicz iterative upstream scheme.
%
%     Number of iterations. This should be in the range 1 - 4. 1 is
%     standard upstream differencing; 3 adds 50% CPU time to POM:
nitera=2;

%     Smoothing parameter. This should preferably be 1, but 0 < sw < 1
%     gives smoother solutions with less overshoot when nitera > 1:
sw=0.50;

%     Index to indicate whether run to start from restart file
%     (nread=0: no restart input file; nread=1: restart input file):
nread=0;

%     External (2-D) time step (secs.) according to CFL:
dte=6.0;

%     <Internal (3-D) time step>/<External (2-D) time step> (dti/dte; dimensionless):
isplit=30;

%     Date and time of start of initial run of model in format (i.e. UDUNITS convention)
%
%       YYYY-MM-DD HH:MM:SS <+/->HH:MM
%
%     where "<+/->HH:MM" is the time zone (positive eastwards from
%     Coordinated Universal Time). NOTE that the climatological time
%     axis (i.e. beginning of year zero, which does not exist in the
%     real-world calendar) has been used here. Insert your own date
%     and time as required:
time_start='2000-01-01 00:00:00 +00:00';


days=0.0625;       % run duration in days
% days=1.5;       % run duration in days

prtd1=0.0125;     % Initial print interval (days)

prtd2=1.0;         % Final print interval (days)

swtch=1000.0;      % Time to switch from prtd1 to prtd2

iskp=4;             % Printout skip interval in i

jskp=3;             % Printout skip interval in j

%     Logical for inertial ramp (true if inertial ramp to be applied
%     to wind stress and baroclinic forcing, otherwise false)
lramp=false;

%     Reference density (recommended values: 1025 for seawater, 1000 for freswater; S.I. units):
rhoref=1025.0;

tbias=0.0;         % Temperature bias (deg. C)

sbias=0.0;         % Salinity bias

grav=9.8060;       % gravity constant (S.I. units)

kappa=0.40;        % von Karman's constant

z0b=.010;          % Bottom roughness (metres)

cbcmin=.00250;     % Minimum bottom friction coeff.

cbcmax=1.0;        % Maximum bottom friction coeff.

horcon=0.20;       % Smagorinsky diffusivity coeff.

%     Inverse horizontal turbulent Prandtl number (ah/am; dimensionless):
%
%     NOTE that tprni=0.0 yields zero horizontal diffusivity!
tprni=.20;

%     Background viscosity used in subroutines profq, proft, profu and profv (S.I. units):
umol=2.e-5;

%     Maximum depth used in radiation boundary condition in subroutine bcond (metres):
hmax=4500.0;

%     Maximum magnitude of vaf (used in check that essentially tests for CFL violation):
vmaxl=100.0;

%     Maximum allowable value of:
%
%       <difference of depths>/<sum of depths>
%
%     for two adjacent cells (dimensionless). This is used in subroutine
%     slpmax. If >= 1, then slpmax is not applied:
slmax=2.0;


%     Water type, used in subroutine proft.
%       ntp    Jerlov water type
%        1            i
%        2            ia
%        3            ib
%        4            ii
%        5            iii
ntp=2;

%     Surface temperature boundary condition, used in subroutine proft:
%       nbct   prescribed    prescribed   short wave
%              temperature      flux      penetration
%        1        no           yes           no
%        2        no           yes           yes
%        3        yes          no            no
%        4        yes          no            yes
%
nbct=1;

%     Surface salinity boundary condition, used in subroutine proft:
%       nbcs   prescribed    prescribed
%               salinity      flux
%        1        no           yes
%        3        yes          no
%     NOTE that only 1 and 3 are allowed for salinity.
nbcs=1;

%     Step interval during which external (2-D) mode advective terms are
%     not updated (dimensionless):
ispadv=5;

%     Constant in temporal filter used to prevent solution splitting
%     (dimensionless):
smoth=0.100;

%     Weight used for surface slope term in external (2-D) dynamic
%     equation (a value of alpha = 0.0 is perfectly acceptable, but the
%     value, alpha=.2250 permits a longer time step):
alpha=0.2250;

%     Initial value of aam:
aam_init=500.0;

% ramp is not assigned in fortran pom2k before used in barop
ramp=0.0;

tatm = 0.0;
satm = 0.0;
io = zeros(100);
jo = zeros(100);
ko = zeros(100);

% gridtype is based on Arakawa grid system;
gridtype = 'C';

%surface forcing scheme 
nsbdy=1;     
%   surfforcing_flag      type       
%         1             seamount     
%         2             box 
%         3             reading from file with constant value
%         4             reading from file with time-varing

npg=1;
%  npg      pressure gradient scheme
%  1       Second order scheme, as originally provide in POM
%  2       Fourth order scheme using the McCalpin method (Berntsen and Oey, Ocean Dynamics, 2010)

% Problem name:
problem='seamount';

% Path of reading data,such as initial value,time-dependent surface forcing data
in_path='./../pre/in/';

%Flags of initial & boundary conditions:
fclim_flag = 0;         % setting initial sclim and tclime(1) or not(0);
wind_flag  = 0;         % reading time-dependent wind force(1) or not(0);
heat_flag  = 0;         % reading time-dependent heat(1) or not(0)
water_flag = 0;         % reading time-dependent salinity (1) or not(0)

bc_flag =0;             % reading time-dependent lateral (1) or not(0)
%     End of input of constants
%***********************************************************************


%     Calculate some constants:
dti=dte*isplit;

dte2=dte*2;
dti2=dti*2;
iend=max(floor(days*24.0*3600.0/dti+0.5),2);

iprint=floor(prtd1*24.0*3600.0/dti+0.5);
iswtch=floor(swtch*24.0*3600.0/dti+0.5);

ispi=1.0/isplit;
isp2i=1.0/(2.0*isplit);
OP = Operator();