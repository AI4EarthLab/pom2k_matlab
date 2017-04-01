function [ua,va,el,et,etf,d,dt,w,d_3d,dt_3d,l,q2b,q2lb,kh,km,kq,aam,q2,q2l,t,s,u,v,rho,drhox,drhoy,drx2d,dry2d] ...
          = update_initial(uab,vab,elb,etb,h,q2b,small,aam,aam_init,w,vfluxf,tb,sb,ub,vb,rmean,ramp)
global kb dz_3d fsm_3d h_3d gs;
ua=uab;     va=vab;     el=elb;     et=etb;     etf=et;
d=h + el;   dt=h + et;


d_3d        = repmat(d,1,1,kb);         dt_3d       =repmat(dt,1,1,kb);
l           = Field(0.1*double(dt_3d),gs,7);    
q2b(:,:,:)  = small;
q2lb        = l.*q2b;                   kh          = l.*sqrt(q2b);
km          = kh;                       kq          = kh;
aam(:,:,:)  = aam_init;                 w(:,:,1)    =vfluxf;

q2          =q2b;                       q2l         =q2lb;
t           =tb;                        s           =sb;
u           =ub;                        v           =vb;

[rho]=dens(s,t,h_3d,fsm_3d);          
[drhox,drhoy] = baropg(rho, rmean, dt_3d, ramp);

drx2d=sum(drhox.*dz_3d, 3);
dry2d=sum(drhoy.*dz_3d, 3);

end