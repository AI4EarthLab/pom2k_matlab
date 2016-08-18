clear all;

im=7;
jm=5;
kb=5;
% im=65;
% jm=49;
% kb=21;
imm1=im-1;
imm2=im-2;
jmm1=jm-1;
jmm2=jm-2;
kbm1=kb-1;
kbm2=kb-2;

alpha          =0.0;dte            =0.0;dti            =0.0;dti2           =0.0;
grav           =0.0;hmax           =0.0;kappa          =0.0;pi             =0.0;
ramp           =0.0;rfe            =0.0;rfn            =0.0;rfs            =0.0;
rfw            =0.0;rhoref         =0.0;sbias          =0.0;slmax          =0.0;
small          =0.0;tbias          =0.0;time           =0.0;tprni          =0.0;
umol           =0.0;vmaxl 		    =0.0;

iint           =0;iprint         =0;iskp           =0;jskp         =0;
kl1            =0;kl2            =0;mode           =0;ntp			=0;

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

%     <Internal (3-D) time step>/<External (2-D) time step>
%     (dti/dte; dimensionless):
isplit=30;

%     Date and time of start of initial run of model in format (i.e.
%     UDUNITS convention)
%
%       YYYY-MM-DD HH:MM:SS <+/->HH:MM
%
%     where "<+/->HH:MM" is the time zone (positive eastwards from
%     Coordinated Universal Time). NOTE that the climatological time
%     axis (i.e. beginning of year zero, which does not exist in the
%     real-world calendar) has been used here. Insert your own date
%     and time as required:
time_start='2000-01-01 00:00:00 +00:00';

days=0.025;       % run duration in days

prtd1=0.0125;     % Initial print interval (days)

prtd2=1.0;         % Final print interval (days)

swtch=1000.0;      % Time to switch from prtd1 to prtd2

iskp=4;             % Printout skip interval in i

jskp=3;             % Printout skip interval in j

%     Logical for inertial ramp (true if inertial ramp to be applied
%     to wind stress and baroclinic forcing, otherwise false)
lramp=false;


%     Reference density (recommended values: 1025 for seawater,
%     1000 for freswater; S.I. units):
rhoref=1025.0;

tbias=0.0;         % Temperature bias (deg. C)

sbias=0.0;         % Salinity bias

grav=9.8060;       % gravity constant (S.I. units)

kappa=0.40;        % von Karman's constant

z0b=.010;          % Bottom roughness (metres)

cbcmin=.00250;     % Minimum bottom friction coeff.

cbcmax=1.0;        % Maximum bottom friction coeff.

horcon=0.20;       % Smagorinsky diffusivity coeff.

%     Inverse horizontal turbulent Prandtl number
%     (ah/am; dimensionless):
%
%     NOTE that tprni=0.0 yields zero horizontal diffusivity!
tprni=.20;

%     Background viscosity used in subroutines profq, proft, profu and
%     profv (S.I. units):
umol=2.e-5;

%     Maximum depth used in radiation boundary condition in subroutine
%     bcond (metres):
hmax=4500.0;

%     Maximum magnitude of vaf (used in check that essentially tests
%     for CFL violation):
vmaxl=100.0;

%     Maximum allowable value of:
%
%       <difference of depths>/<sum of depths>
%
%     for two adjacent cells (dimensionless). This is used in subroutine
%     slpmax. If >= 1, then slpmax is not applied:
slmax=2.0;


%     Integers defining the number of logarithmic layers at the
%     surface and bottom (used by subroutine depth). The number of
%     logarithmic layers are kl1-2 at the surface and kb-kl2-1
%     at the bottom. For no log portions, set kl1=2 and kl2=kb-1:
kl1=6;
kl2=kb-2;

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


dz=zeros(1,kb)       ;dzz=zeros(1,kb)      ;z=zeros(1,kb)        ;zz=zeros(1,kb);

aam2d=zeros(im,jm)   ;advua=zeros(im,jm)   ;advva=zeros(im,jm)   ;adx2d=zeros(im,jm)   ;
ady2d=zeros(im,jm)   ;art=zeros(im,jm)     ;aru=zeros(im,jm)     ;arv=zeros(im,jm)     ;
cbc=zeros(im,jm)     ;cor=zeros(im,jm)     ;d=zeros(im,jm)       ;drx2d=zeros(im,jm)   ;
dry2d=zeros(im,jm)   ;dt=zeros(im,jm)      ;dum=zeros(im,jm)     ;dvm=zeros(im,jm)     ;
dx=zeros(im,jm)      ;dy=zeros(im,jm)      ;east_c=zeros(im,jm)  ;east_e=zeros(im,jm)  ;
east_u=zeros(im,jm)  ;east_v=zeros(im,jm)  ;e_atmos=zeros(im,jm) ;egb=zeros(im,jm)     ;
egf=zeros(im,jm)     ;el=zeros(im,jm)      ;elb=zeros(im,jm)     ;elf=zeros(im,jm)     ;
et=zeros(im,jm)      ;etb=zeros(im,jm)     ;etf=zeros(im,jm)     ;fluxua=zeros(im,jm)  ;
fluxva=zeros(im,jm)  ;fsm=zeros(im,jm)     ;h=zeros(im,jm)       ;north_c=zeros(im,jm) ;
north_e=zeros(im,jm) ;north_u=zeros(im,jm) ;north_v=zeros(im,jm) ;psi=zeros(im,jm)     ;
rot=zeros(im,jm)     ;ssurf=zeros(im,jm)   ;swrad=zeros(im,jm)   ;vfluxb=zeros(im,jm)  ;
tps=zeros(im,jm)     ;tsurf=zeros(im,jm)   ;ua=zeros(im,jm)      ;vfluxf=zeros(im,jm)  ;
uab=zeros(im,jm)     ;uaf=zeros(im,jm)     ;utb=zeros(im,jm)     ;utf=zeros(im,jm)     ;
va=zeros(im,jm)      ;vab=zeros(im,jm)     ;vaf=zeros(im,jm)     ;
vtb=zeros(im,jm)     ;vtf=zeros(im,jm)     ;wssurf=zeros(im,jm)  ;wtsurf=zeros(im,jm)  ;
wubot=zeros(im,jm)   ;wusurf=zeros(im,jm)  ;wvbot=zeros(im,jm)   ;wvsurf=zeros(im,jm) ;


aam=zeros(im,jm,kb)  ;advx=zeros(im,jm,kb) ;advy=zeros(im,jm,kb) ;a=zeros(im,jm,kb)    ;
c=zeros(im,jm,kb)    ;drhox=zeros(im,jm,kb);drhoy=zeros(im,jm,kb);dtef=zeros(im,jm,kb) ;
ee=zeros(im,jm,kb)   ;gg=zeros(im,jm,kb)   ;kh=zeros(im,jm,kb)   ;km=zeros(im,jm,kb)   ;
kq=zeros(im,jm,kb)   ;l=zeros(im,jm,kb)    ;q2b=zeros(im,jm,kb)  ;q2=zeros(im,jm,kb)   ;
q2lb=zeros(im,jm,kb) ;q2l=zeros(im,jm,kb)  ;rho=zeros(im,jm,kb)  ;rmean=zeros(im,jm,kb);
sb=zeros(im,jm,kb)   ;sclim=zeros(im,jm,kb);s=zeros(im,jm,kb)    ;tb=zeros(im,jm,kb)   ;
tclim=zeros(im,jm,kb);t=zeros(im,jm,kb)    ;ub=zeros(im,jm,kb)   ;uf=zeros(im,jm,kb)   ;
u=zeros(im,jm,kb)    ;vb=zeros(im,jm,kb)   ;vf=zeros(im,jm,kb)   ;v=zeros(im,jm,kb)    ;
w=zeros(im,jm,kb)    ;zflux=zeros(im,jm,kb);

ele=zeros(jm)        ;eln=zeros(im)        ;els=zeros(im)        ;elw=zeros(jm)        ;
sbe=zeros(jm,kb)     ;sbn=zeros(im,kb)     ;sbs=zeros(im,kb)     ;sbw=zeros(jm,kb)     ;
tbe=zeros(jm,kb)     ;tbn=zeros(im,kb)     ;tbs=zeros(im,kb)     ;tbw=zeros(jm,kb)     ;
uabe=zeros(jm)       ;uabw=zeros(jm)       ;ube=zeros(jm,kb)     ;ubw=zeros(jm,kb)     ;
vabn=zeros(im)       ;vabs=zeros(im)       ;vbn=zeros(im,kb)     ;vbs=zeros(im,kb);



if(iproblem ~= 3)
    [kdz,z,zz,dz,dzz]=depth(z,zz,dz,dzz,kb,kl1,kl2);
end



if(iproblem == 1)
    
    [dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs] = seamount(dx,dy,cor,...
        east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
        h,art,aru,arv,fsm,dum,dvm,...
        tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
        aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
        uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,...
        e_atmos,aam,im,jm,kb,imm1,jmm1,kbm1,slmax,zz,tbias,sbias,grav,rhoref);

    
elseif(iproblem == 2)
   [dx,dy,cor,...
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
        e_atmos,aam,im,jm,kb,imm1,jmm1,kbm1,slmax,zz,tbias,sbias,grav,rhoref,rot,tatm,satm,vfluxf);
     
elseif(iproblem ==3)
    [dx,dy,cor,...
    east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
    h,art,aru,arv,fsm,dum,dvm,...
    tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
    aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
    uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,rot,els,eln,vabs,vabn,ubw,ube,vbs,vbn] = file2ic(dx,dy,cor,...
    east_c,north_c,east_e,north_e,east_u,north_u,east_v,north_v,...
    h,art,aru,arv,fsm,dum,dvm,...
    tb,sb,tclim,sclim,ub,uab,elb,etb,dt,...
    aam2d,rho,rmean,rfe,rfw,rfn,rfs,...
    uabw,uabe,ele,elw,tbe,tbw,sbe,sbw,tbn,tbs,sbn,sbs,...
    e_atmos,aam,im,jm,kb,imm1,jmm1,kbm1,slmax,zz,tbias,sbias,grav,rhoref,rot,els,eln,vabs,vabn,ubw,ube,vbs,vbn);
    
else
    fprintf('Invalid value of iproblem ..... program terminated\n');
    return
end

%     Inertial period for temporal filter:
period=(2.0*pi)/abs(cor(floor(im/2),floor(jm/2)))/86400.0;

%     Initialise time:
time0=0.0;
time=0.0;
%     Initial conditions:
%     NOTE that lateral thermodynamic boundary conditions are often set
%     equal to the initial conditions and are held constant thereafter.
%    Users can of course create variable boundary conditions.

ua=uab;
va=vab;
el=elb;
et=etb;
etf=et;
d=h + el;
dt=h + et;
w(:,:,1)=vfluxf;

% % % for k=1:kb
% % %     for j=1:jm
% % %         for i=1:im
% % %             l(i,j,k)=0.1*dt(i,j);
% % %             q2b(i,j,k)=small;
% % %             q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k);
% % %             kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k));
% % %             km(i,j,k)=kh(i,j,k);
% % %             kq(i,j,k)=kh(i,j,k);
% % %             aam(i,j,k)=aam_init;
% % %         end
% % %     end
% % % end

for k=1:kb
    l(:,:,k)    = 0.1*dt(:,:);
    q2b(:,:,k)  = small;
    q2lb(:,:,k) = l(:,:,k).*q2b(:,:,k);
    kh(:,:,k)   = l(:,:,k).*sqrt(q2b(:,:,k));
    km(:,:,k)   = kh(:,:,k);
    kq(:,:,k)   = kh(:,:,k);
    aam(:,:,k)  = aam_init;
end

% % % for k=1:kbm1
% % %     for i=1:im
% % %         for j=1:jm
% % %             q2(i,j,k)=q2b(i,j,k);
% % %             q2l(i,j,k)=q2lb(i,j,k);
% % %             t(i,j,k)=tb(i,j,k);
% % %             s(i,j,k)=sb(i,j,k);
% % %             u(i,j,k)=ub(i,j,k);
% % %             v(i,j,k)=vb(i,j,k);
% % %         end
% % %     end
% % % end

for k=1:kbm1
    q2(:,:,k) =q2b(:,:,k);
    q2l(:,:,k)=q2lb(:,:,k);
    t(:,:,k)  =tb(:,:,k);
    s(:,:,k)  =sb(:,:,k);
    u(:,:,k)  =ub(:,:,k);
    v(:,:,k)  =vb(:,:,k);
end

[rho]=dens(s,t,rho,im,jm,kbm1,tbias,sbias,grav,rhoref,zz,h,fsm);

[rho,drhox,drhoy] = baropg(rho,drhox,drhoy,...
    im,jm,imm1,jmm1,kb,kbm1,grav,...
    zz,dt,dum,dvm,ramp,rmean,dx,dy);



for k=1:kbm1
    for j=1:jm
        for i=1:im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k);
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k);
        end
    end
end
%
%     Calculate bottom friction coefficient:
%
for j=1:jm
    for i=1:im
        cbc(i,j)=(kappa/log((1.e0+zz(kbm1))*h(i,j)/z0b))^2;
        cbc(i,j)=max(cbcmin,cbc(i,j));
        %
        %     If the following is invoked, then it is probable that the wrong
        %     choice of z0b or vertical spacing has been made:
        %
        cbc(i,j)=min(cbcmax,cbc(i,j));
    end
end

%
%     Calculate external (2-D) CFL time step:
%
for j=1:jm
    for i=1:im
        tps(i,j)=0.5e0/sqrt(1.e0/dx(i,j)^2+1.e0/dy(i,j)^2)...
            /sqrt(grav*(h(i,j)+small))*fsm(i,j);
    end
end
d = h+el;
dt = h+et;
time = time0;

%==========================================
%           begin internal (3-D) mode
%
%==========================================

for  iint=1:iend
    time=dti*iint*1.0/86400+time0;
    if(lramp~=0)
        ramp = time/period;
        if(ramp>1.0)
            ramp=1.0;
        end
    else
        ramp=1.0;
    end
    %
    %       write(6,2) mode,iint,time
    %   2   format(' mode,iint,time =',2i5,f9.2)
    %
    %-----------------------------------------------------------------------
    %
    %     Set time dependent, surface and lateral boundary conditions.
    %     The latter will be used in subroutine bcond. Users may
    %     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
    %     wssurf, swrad and vflux.
    %
    %     Introduce simple wind stress. Value is negative for westerly or
    %     southerly winds. The following wind stress has been tapered
    %     along the boundary to suppress numerically induced oscilations
    %     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
    %     To make a healthy surface Ekman layer, it would be well to set
    %     kl1=9.
    %
    for j=2:jmm1
        for i=2:imm1
            if(iproblem~=3) % constant wind read in file2ic
                
    %           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1));
                wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))  ...
                    *0.25*(dvm(i,j+1)+dvm(i-1,j+1)     ...
                    +dvm(i-1,j)+dvm(i,j));
                % --- no wind ----
                %           wusurf(i,j)=0.e0;
                wvsurf(i,j)=0.e0;
            end
            
            e_atmos(i,j)=0.e0;
            vfluxf(i,j)=0.e0;
            %
            %     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
            %     the sea surface. See calculation of elf(i,j) below and subroutines
            %     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
            %     is no net flow across lateral boundaries, the basin volume will be
            %     constant; if also vflux(i,j).ne.0, then, for example, the average
            %     salinity will change and, unrealistically, so will total salt.
            %
            w(i,j,1)=vfluxf(i,j);
            %
            %     Set wtsurf to the sensible heat, the latent heat (which involves
            %     only the evaporative component of vflux) and the long wave
            %     radiation:
            %
            wtsurf(i,j)=0.0;
            %
            %     Set swrad to the short wave radiation:
            %
            swrad(i,j)=0.0;
            %
            %     To account for change in temperature of flow crossing the sea
            %     surface (generally quite small compared to latent heat effect)
            %
            tatm=t(i,j,1)+tbias;    % an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias);
            %
            %     Set the salinity of water vapor/precipitation which enters/leaves
            %     the atmosphere (or e.g., an ice cover)
            %
            satm=0.0  ;
            wssurf(i,j)=  vfluxf(i,j)*(satm-s(i,j,1)-sbias)  ;
            %
        end
    end
    %
    %-----------------------------------------------------------------------
    %
    %     Set lateral viscosity:
    %
    %     If mode=2 then initial values of aam2d are used. If one wishes
    %     to use Smagorinsky lateral viscosity and diffusion for an
    %     external (2-D) mode calculation, then appropiate code can be
    %     adapted from that below and installed just before the end of the
    %     "if(mode.eq.2)" loop in subroutine advave.
    %
    %     %alculate Smagorinsky lateral viscosity:
    %
    %       ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
    %                                     +.5*(du/dy+dv/dx)**2) )
    %
    if(mode~=2)
        [a,c,ee,advx,advy]=...
            advct(a,c,ee,advx,advy,...
            u,v,dx,dy,dt,aam,ub,vb,aru,arv,im,jm,kb,imm1,jmm1,kbm1);
        
        [rho,drhox,drhoy] = baropg(rho,drhox,drhoy,...
            im,jm,imm1,jmm1,kb,kbm1,grav,...
            zz,dt,dum,dvm,ramp,rmean,dx,dy);
        
        
        for k=1:kbm1
            for j=2:jmm1
                for i=2:imm1
                    aam(i,j,k)=horcon*dx(i,j)*dy(i,j) ...
                        *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))^2     ...
                        +((v(i,j+1,k)-v(i,j,k))/dy(i,j))^2     ...
                        +0.5*(0.25*(u(i,j+1,k)+u(i+1,j+1,k)     ...
                        -u(i,j-1,k)-u(i+1,j-1,k))     ...
                        /dy(i,j)     ...
                        +0.25*(v(i+1,j,k)+v(i+1,j+1,k)     ...
                        -v(i-1,j,k)-v(i-1,j+1,k))     ...
                        /dx(i,j))^2);
                end
            end
        end
        %
        %     Form vertical averages of 3-D fields for use in external (2-D)
        %     mode:
        %
        
        adx2d=zeros(im,jm);
        ady2d=zeros(im,jm);
        drx2d=zeros(im,jm);
        dry2d=zeros(im,jm);
        aam2d=zeros(im,jm);
        
        %
        for k=1:kbm1
            for j=1:jm
                for i=1:im
                    adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k);
                    ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k);
                    drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k);
                    dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k);
                    aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k);
                end
            end
        end
        
        %
        [tps,advua,advva,fluxua,fluxva,wubot,wvbot,tps] = advave(tps,advua,advva,fluxua,fluxva,wubot,wvbot,tps,...
            mode,im,jm,imm1,jmm1,aam2d,...
            uab,vab,dx,dy,ua,va,cbc,aru,arv,d);
        
        %
        for j=1:jm
            for i=1:im
                adx2d(i,j)=adx2d(i,j)-advua(i,j);
                ady2d(i,j)=ady2d(i,j)-advva(i,j);
            end
        end
        %
    end
    %
    for j=1:jm
        for i=1:im
            egf(i,j)=el(i,j)*ispi;
        end
    end
    %
    for j=1:jm
        for i=2:im
            utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i;
        end
    end
    for j=2:jm
        for i=1:im
            vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i;
        end
    end
    
    
    %
    %-----------------------------------------------------------------------
    
    
    
    
    
    for iext=1:isplit    % Begin external (2-D) mode
        %
        %         write(6,3) iext,time
        %   3     format(' iext,time =',i5,f9.2)
        %
        for j=2:jm
            for i=2:im
                fluxua(i,j)=0.25*(d(i,j)+d(i-1,j))  ...
                    *(dy(i,j)+dy(i-1,j))*ua(i,j);
                fluxva(i,j)=0.25*(d(i,j)+d(i,j-1))     ...
                    *(dx(i,j)+dx(i,j-1))*va(i,j);
            end
        end
        %
        %     NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
        %     with pom98.f. See also modifications to subroutine vertvl.
        %
        
        
        for j=2:jmm1
            for i=2:imm1
                elf(i,j)=elb(i,j)     ...
                    +dte2*(-(fluxua(i+1,j)-fluxua(i,j)     ...
                    +fluxva(i,j+1)-fluxva(i,j))/art(i,j)     ...
                    -vfluxf(i,j));
            end
        end
        %
        
        
        [elf,uaf,vaf,uf,vf,w] = bcond(1,elf,uaf,vaf,uf,vf,w,...
            im,jm,kb,imm1,jmm1,kbm1,...
            fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
            dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        
        
        
        if(mod(iext,ispadv)==0)
            [tps,advua,advva,fluxua,fluxva,wubot,wvbot,tps] = advave(tps,advua,advva,fluxua,fluxva,wubot,wvbot,tps,...
                mode,im,jm,imm1,jmm1,aam2d,...
                uab,vab,dx,dy,ua,va,cbc,aru,arv,d);
          
        end
        
        
        
        
        %
        for j=2:jmm1
            for i=2:im
                uaf(i,j)=adx2d(i,j)+advua(i,j)     ...
                    -aru(i,j)*0.25     ...
                    *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))     ...
                    +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))     ...
                    +0.25*grav*(dy(i,j)+dy(i-1,j))     ...
                    *(d(i,j)+d(i-1,j))     ...
                    *((1.e0-2.0*alpha)     ...
                    *(el(i,j)-el(i-1,j))     ...
                    +alpha*(elb(i,j)-elb(i-1,j)     ...
                    +elf(i,j)-elf(i-1,j))     ...
                    +e_atmos(i,j)-e_atmos(i-1,j))     ...
                    +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j));
            end
        end
        %
        for j=2:jmm1
            for i=2:im
                uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))     ...
                    *aru(i,j)*uab(i,j)     ...
                    -4.e0*dte*uaf(i,j))     ...
                    /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))     ...
                    *aru(i,j));
            end
        end
        
        %
        for j=2:jm
            for i=2:imm1
                vaf(i,j)=ady2d(i,j)+advva(i,j)     ...
                    +arv(i,j)*0.25     ...
                    *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))     ...
                    +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))     ...
                    +0.25*grav*(dx(i,j)+dx(i,j-1))     ...
                    *(d(i,j)+d(i,j-1))     ...
                    *((1.e0-2.0*alpha)*(el(i,j)-el(i,j-1))     ...
                    +alpha*(elb(i,j)-elb(i,j-1)     ...
                    +elf(i,j)-elf(i,j-1))     ...
                    +e_atmos(i,j)-e_atmos(i,j-1))     ...
                    +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j));
            end
        end
        %
        for j=2:jm
            for i=2:imm1
                vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))     ...
                    *vab(i,j)*arv(i,j)     ...
                    -4.e0*dte*vaf(i,j))     ...
                    /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))     ...
                    *arv(i,j));
            end
        end
        %
        
        [elf,uaf,vaf,uf,vf,w] = bcond(2,elf,uaf,vaf,uf,vf,w,...
            im,jm,kb,imm1,jmm1,kbm1,...
            fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
            dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
        %
        if(iext==(isplit-2))
            for j=1:jm
                for i=1:im
                    etf(i,j)=0.25*smoth*elf(i,j);
                end
            end
            %
        elseif(iext==(isplit-1))
            %
            for j=1:jm
                for i=1:im
                    etf(i,j)=etf(i,j)+0.5*(1.-.5e0*smoth)*elf(i,j);
                end
            end
            %
        elseif(iext==isplit)
            %
            for j=1:jm
                for i=1:im
                    etf(i,j)=(etf(i,j)+0.5*elf(i,j))*fsm(i,j);
                end
            end
            %
        end
        %
        %     Stop if velocity condition violated (generally due to %FL
        %     criterion not being satisfied):
        %
        
        vamax=0.0;
        %
        for j=1:jm
            for i=1:im
                if(abs(vaf(i,j))>=vamax)
                    vamax=abs(vaf(i,j));
                    imax=i;
                    jmax=j;
                end
            end
        end
        %
        if(vamax<=vmaxl)
            %
            %     Apply filter to remove time split and reset time sequence:
            %
            for j=1:jm
                for i=1:im
                    ua(i,j)=ua(i,j)     ...
                        +0.5*smoth*(uab(i,j)-2.0*ua(i,j)+uaf(i,j));
                    va(i,j)=va(i,j)     ...
                        +0.5*smoth*(vab(i,j)-2.0*va(i,j)+vaf(i,j));
                    el(i,j)=el(i,j)     ...
                        +0.5*smoth*(elb(i,j)-2.0*el(i,j)+elf(i,j));
                    elb(i,j)=el(i,j);
                    el(i,j)=elf(i,j);
                    d(i,j)=h(i,j)+el(i,j);
                    uab(i,j)=ua(i,j);
                    ua(i,j)=uaf(i,j);
                    vab(i,j)=va(i,j);
                    va(i,j)=vaf(i,j);
                end
            end
            %
            if(iext~=isplit)
                for j=1:jm
                    for i=1:im
                        egf(i,j)=egf(i,j)+el(i,j)*ispi;
                    end
                end
                for j=1:jm
                    for i=2:im
                        utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i;
                    end
                end
                for j=2:jm
                    for i=1:im
                        vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i;
                    end
                end
            end
            %
        end
        %
    end %
    
    
    
    %===========================================
    %End of external (2-D) mode
    %
    %
    %=============================================
    
    if(vamax<=vmaxl)
        
        %
        %     continue with internal (3-D) mode calculation:
        %
        if((iint~= 1|| time0~=0.e0)&&mode~=2)
            %
            %     Adjust u(z) and v(z) such that depth average of (u,v) = (ua,va):
            %
            tps=zeros(im,jm);
            %
            for k=1:kbm1
                for j=1:jm
                    for i=1:im
                        
                        tps(i,j)=tps(i,j)+u(i,j,k)*dz(k);
                    end
                end
            end
            
            %
            for k=1:kbm1
                for j=1:jm
                    for i=2:im
                        u(i,j,k)=(u(i,j,k)-tps(i,j))+     ...
                            (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j));
                    end
                end
            end
            
            
            %
            
            tps =zeros(im,jm);
            %
            for k=1:kbm1
                for j=1:jm
                    for i=1:im
                        tps(i,j)=tps(i,j)+v(i,j,k)*dz(k);
                    end
                end
            end
            
            for k=1:kbm1
                for j=2:jm
                    for i=1:im
                        v(i,j,k)=(v(i,j,k)-tps(i,j))+ ...
                            (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1));
                    end
                end
            end
            
            
            %     vertvl calculates w from u, v, dt (h+et), etf and etb:
            %
            
            [a,c,...
                w]=vertvl(a,c,...
                w,dx,dy,dz,dt,u,v,vfluxb,vfluxf,etf,etb,dti2,im,jm,imm1,jmm1,kbm1);
            
            [elf,uaf,vaf,uf,vf,w] = bcond(5,elf,uaf,vaf,uf,vf,w,...
                im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
            
            
            
            vf=zeros(im,jm,kb);
            uf=zeros(im,jm,kb);
            %            %
            %            %     %alculate q2f and q2lf using uf, vf, a and c as temporary
            %            %     variables:
            %            %
            %
            
            [q2b,q2,uf,a,c]=advq(q2b,q2,uf,a,c,...
                dt,dx,dy,dz,u,v,w,aam,h,dum,dvm,art,etb,etf,im,jm,imm1,jmm1,kbm1,dti2);
            
            
            [q2lb,q2l,vf,a,c]=advq(q2lb,q2l,vf,a,c,...
                dt,dx,dy,dz,u,v,w,aam,h,dum,dvm,art,etb,etf,im,jm,imm1,jmm1,kbm1,dti2);
            
            
            
            
            [a,c,tps,dtef,...
                ee,gg,l,kq,km,kh,...
                uf,vf,q2b,q2lb,a,c]=profq(a,c,tps,dtef,....
                ee,gg,l,kq,km,kh,...
                uf,vf,q2,q2b,q2lb,a,c,...
                h,etf,dti2,umol,dzz,grav,rho,kappa,u,v,dt,small,fsm,im,jm,kb,imm1,jmm1,kbm1,tbias,sbias,dz,...
                wusurf,wubot,wvsurf,wvbot,t,s,rhoref,zz,z);
            
            
            [elf,uaf,vaf,uf,vf,w] = bcond(6,elf,uaf,vaf,uf,vf,w,...
                im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
            
            for k=1:kb
                for j=1:jm
                    for i=1:im
                        q2(i,j,k)=q2(i,j,k)    ...
                            +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)     ...
                            -2.e0*q2(i,j,k));
                        q2l(i,j,k)=q2l(i,j,k)     ...
                            +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)     ...
                            -2.e0*q2l(i,j,k));
                        q2b(i,j,k)=q2(i,j,k);
                        q2(i,j,k)=uf(i,j,k);
                        q2lb(i,j,k)=q2l(i,j,k);
                        q2l(i,j,k)=vf(i,j,k);
                    end
                end
            end
            %
            %      calculate tf and sf using uf, vf, a and c as temporary variables:
            %
            if(mode~=4)
                %
                if(nadv==1)
                    
                    [tb,t,tclim,uf,a,c,zflux]=advt1(tb,t,tclim,uf,a,c,...
                        zflux,dt,u,v,aam,tprni,...
                        dum,dvm,w,art,etb,etf,dti2,dx,dy,dz,h,im,jm,kb,imm1,jmm1,kbm1);
                    
                    [sb,s,sclim,vf,a,c,zflux]=advt1(sb,s,sclim,vf,a,c,...
                        zflux,dt,u,v,aam,tprni,...
                        dum,dvm,w,art,etb,etf,dti2,dx,dy,dz,h,im,jm,kb,imm1,jmm1,kbm1);
                    
                elseif(nadv==2)
                    %
                   
                    [tb,t,tclim,uf,a,c,nitera,sw,...
                     zflux] = advt2(tb,t,tclim,uf,a,c,nitera,sw,...
                                                 zflux,...
                                                 im,jm,kb,imm1,jmm1,kbm1,dti2,...
                                                 etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);
                    [sb,s,sclim,vf,a,c,nitera,sw,...
                     zflux] = advt2(sb,s,sclim,vf,a,c,nitera,sw,...
                                                 zflux,...
                                                 im,jm,kb,imm1,jmm1,kbm1,dti2,...
                                                 etb,etf,w,art,fsm,dt,aam,tprni,h,dum,dvm,dx,dy,u,v,aru,arv,dz,dzz);      
                    %
                else
                    %
                    
                    disp 'Invalid value for nadv ..... '
                    disp 'program terminated'
                    
                    return
                    %
                end
                %
                
                [uf,wtsurf,tsurf,nbct,tps,...
                    a,c,ee,gg] = proft(uf,wtsurf,tsurf,nbct,tps,...
                    a,c,ee,gg,...
                    h,etf,dti2,dz,dzz,swrad,ntp,im,jm,kb,kbm1,kbm2,kh,umol);
                
                %
                [vf,wssurf,ssurf,nbcs,tps,...
                    a,c,ee,gg] = proft(vf,wssurf,ssurf,nbcs,tps,...
                    a,c,ee,gg,...
                    h,etf,dti2,dz,dzz,swrad,ntp,im,jm,kb,kbm1,kbm2,kh,umol);
                %
                [elf,uaf,vaf,uf,vf,w] = bcond(4,elf,uaf,vaf,uf,vf,w,...
                    im,jm,kb,imm1,jmm1,kbm1,...
                    fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                    dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
                for k=1:kb
                    for j=1:jm
                        for i=1:im
                            t(i,j,k)=t(i,j,k)     ...
                                +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)     ...
                                -2.e0*t(i,j,k));
                            s(i,j,k)=s(i,j,k)     ...
                                +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)     ...
                                -2.e0*s(i,j,k));
                            tb(i,j,k)=t(i,j,k);
                            t(i,j,k)=uf(i,j,k);
                            sb(i,j,k)=s(i,j,k);
                            s(i,j,k)=vf(i,j,k);
                        end
                    end
                end
                %
                %
                %
                [rho]=dens(s,t,rho,...
                    im,jm,kbm1,tbias,sbias,grav,rhoref,zz,h,fsm);
                
                %
            end  % end if
            %
            %     calculate uf and vf:
            %
            [uf] = advu(advx,aru,dz,cor,dt,e_atmos,dy,drhox,h,dti2,ub,uf,u,v,w,im,jm,kb,imm1,jmm1,kbm1,grav,egf,egb,etf,etb);
            [vf] = advv(advy,arv,dz,cor,dt,e_atmos,dx,drhoy,h,dti2,vb,vf,u,v,w,im,jm,kb,imm1,jmm1,kbm1,grav,egf,egb,etf,etb);
            [a,c,ee,gg,tps,uf,wubot] = profu(a,c,ee,gg,tps,uf,wubot,...
                etf,h,km,dti2,umol,dz,dzz,wusurf,cbc,dum,im,jm,kb,imm1,jmm1,kbm1,kbm2,ub,vb);
            
            [a,c,ee,gg,tps,vf,wvbot] = profv(a,c,ee,gg,tps,vf,wvbot,...
                dvm,dz,dzz,im,jm,kb,imm1,jmm1,kbm1,kbm2,...
                km,cbc,ub,vb,umol,wvsurf,h,etf,dti2) ;
            
            [elf,uaf,vaf,uf,vf,w] = bcond(3,elf,uaf,vaf,uf,vf,w,...
                im,jm,kb,imm1,jmm1,kbm1,...
                fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz);
            %
            
            tps = zeros(im,jm);
            %
            for k=1:kbm1
                for j=1:jm
                    for i=1:im
                        tps(i,j)=tps(i,j)     ...
                            +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k);
                    end
                end
            end
            %
            for k=1:kbm1
                for j=1:jm
                    for i=1:im
                        u(i,j,k)=u(i,j,k)     ...
                            +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)     ...
                            -2.e0*u(i,j,k)-tps(i,j));
                    end
                end
            end
            %
            tps = zeros(im,jm);
            %
            for k=1:kbm1
                for j=1:jm
                    for i=1:im
                        tps(i,j)=tps(i,j)     ...
                            +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k);
                    end
                end
            end
            %
            for k=1:kbm1
                for j=1:jm
                    for i=1:im
                        v(i,j,k)=v(i,j,k)     ...
                            +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)     ...
                            -2.e0*v(i,j,k)-tps(i,j));
                    end
                end
            end
            %
            
            ub=u;
            u=uf;
            vb=v;
            v=vf;
            %
        end  % end if
        %        %
        egb=egf;
        etb=et;
        et=etf;
        dt=h+et;
        utb=utf;
        vtb=vtf;
        vfluxb=vfluxf;
        %
    end   % end if
    
  

    %
    %     Beginning of print section:
    %
    if(iint>=iswtch)
        iprint=floor(prtd2*24.e0*3600.e0/dti + 0.5);
    end
    %
    if(mod(iint,iprint)==0 || vamax>=vmaxl)
        fprintf('**************************************************\n');
        fprintf('**************************************************\n');
        fprintf('time = %.6f, iint = %d, iexit = %d, iprint = %d\n',time,iint,iext,iprint);
        
        %
        %     Select print statements in printall as desired:
        %
        %printall
        %
        vtot=0.e0;
        atot=0.e0;
        taver=0.e0;
        saver=0.e0;
        eaver=0.e0;
        for k=1:kbm1
            for j=1:jm
                for i=1:im
                    darea=dx(i,j)*dy(i,j)*fsm(i,j);
                    dvol=darea*dt(i,j)*dz(k);
                    vtot=vtot+dvol;
                    taver=taver+tb(i,j,k)*dvol;
                    saver=saver+sb(i,j,k)*dvol;
                end
            end
        end
        %
        for j=1:jm
            for i=1:im
                darea=dx(i,j)*dy(i,j)*fsm(i,j);
                atot=atot+darea;
                eaver=eaver+et(i,j)*darea;
            end
        end
        %
        taver=taver/vtot;
        saver=saver/vtot;
        eaver=eaver/atot;
        tsalt=(saver+sbias)*vtot;
      
        fprintf('vtot = %.6f,atot = %.6f\n',vtot,atot);
        fprintf('eaver = %.6f,taver = %.6f,saver=%.6f,saver = %.6f,tsalt = %.6f\n',eaver,taver,saver,tsalt,tsalt);
        
        
        if(vamax>vmaxl)
            
            fprintf('time = %.6f, iint = %.6f, iexit = %.6f, iprint = %.6f\n',time,iint,iext,iprint);
           
            printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s,rho,aam,km,kb,mode,dt,zz,z);
           
            disp        '************************************************'
            disp        '************ abnormal job end ******************'
            disp        '************* user terminated ******************'
            disp        '************************************************'
            
            fprintf('vamax = %d, imax = %d ,jmax = %d \n',vmax,imax,jmax);
            return;
            %
        end
        %
    end
    %
    %     End of print section
    %
    %-----------------------------------------------------------------------
    %
    %  End of internal (3-D) mode
    %
    %-----------------------------------------------------------------------
    %
    
    
end
%internal loop ended
%write(6,4) time,iint,iext,iprint

fprintf('time = %d, iint = %d ,iext = %d , iprint = %d \n',time,iint,iext,iprint);
%
%     Set levels for output:
%
ko(1)=1;
ko(2)=2;
ko(3)=kb/2;
ko(4)=kb-1;
ko(5)=kb;
%
%      prxyz('Vertical velocity, w                    ',
%    ...           time,w       ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
%
%      prxyz('Turbulent kinetic energy x 2, q2        ',
%    ...           time,q2      ,im,iskp,jm,jskp,kb,ko,5,-1.e0)
%
%     Save this data for a seamless restart:
%

%    write_to_file(im,jm,kb,time,wubot,wvbot,aam2d,ua,uab,va,vab,el,elb,et,etb,egb,...
%    utb,vtb,u,ub,w,v,vb,t,tb,s,sb,rho,adx2d,ady2d,advua,advva,...
%    km,kh,kq,l,q2,q2b,aam,q2l,q2lb);

    printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s,rho,aam,km,kb,mode,dt,zz,z);
    
%
%
%
%     End of main program