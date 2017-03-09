function [period,time,time0,ua,va,el,et,etf,d,dt,w,d_3d,dt_3d,l,q2b,q2lb,kh,km,kq,aam,q2,q2l,t,s,u,v,rho, ...
          drhox,drhoy,drx2d,dry2d,cbc,tps] = init_time_condition(uab,vab,elb,etb,vfluxf,tb,sb,ub,vb,rmean,ramp,q2b,aam,w)
global kb small im jm pi h aam_init h_3d fsm_3d dz_3d zz kappa kbm1 z0b cbcmin cbcmax ...
       dx dy grav fsm cor
period=(2.0*pi)/abs(cor(floor(im/2),floor(jm/2)))/86400.0;

time0=0.0;  time =0.0;
ua=uab;     va=vab;     el=elb;     et=etb;     etf=et;
d=h + el;   dt=h + et;
w(:,:,1)=vfluxf;

d_3d=repmat(d,1,1,kb);      dt_3d=repmat(dt,1,1,kb);
l           = 0.1*dt_3d;    q2b(:,:,:)  = small;
q2lb        = l.*q2b;       kh          = l.*sqrt(q2b);
km          = kh;           kq          = kh;
aam(:,:,:)  = aam_init;

q2 =q2b;                    q2l=q2lb;
t  =tb;                     s  =sb;
u  =ub;                     v  =vb;

[rho]=dens(s,t,h_3d,fsm_3d);          
[drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp);

drx2d=sum(drhox.*dz_3d, 3);
dry2d=sum(drhoy.*dz_3d, 3);

cbc=(kappa./log((1.0+zz(kbm1))*h/z0b)).^2;
cbc=max(cbcmin,cbc);        cbc=min(cbcmax,cbc);

%     Calculate external (2-D) CFL time step:
tps=0.5./sqrt(1.0./dx.^2+1.0./dy.^2) ./ sqrt(grav*(h+small)) .* fsm;
time = time0;
return
