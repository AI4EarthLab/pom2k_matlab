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

d_3d=repmat(d,1,1,kb);
dt_3d=repmat(dt,1,1,kb);


l           = 0.1*dt_3d;
q2b(:,:,:)  = small;
q2lb        = l.*q2b;
kh          = l.*sqrt(q2b);
km          = kh;
kq          = kh;
aam(:,:,:)  = aam_init;

q2 =q2b;
q2l=q2lb;
t  =tb;
s  =sb;
u  =ub;
v  =vb;

[rho]=dens(s,t,h_3d,fsm_3d);

                      
 [drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp);

drx2d=sum(drhox.*dz_3d, 3);
dry2d=sum(drhoy.*dz_3d, 3);

%     Calculate bottom friction coefficient:
cbc=(kappa./log((1.0+zz(kbm1))*h/z0b)).^2;
cbc=max(cbcmin,cbc);
%     If the following is invoked, then it is probable that the wrong
%     choice of z0b or vertical spacing has been made:
cbc=min(cbcmax,cbc);

%     Calculate external (2-D) CFL time step:
tps=0.5./sqrt(1.0./dx.^2+1.0./dy.^2) ./ sqrt(grav*(h+small)) .* fsm;

time = time0;
