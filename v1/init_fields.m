
aam2d   = create_field(aam2d,   gs, 3);

adx2d   = create_field(adx2d,   gs, 2);
ady2d   = create_field(ady2d,   gs, 1);
art_3d  = create_field(art_3d,  gs, 3);
aru_3d  = create_field(aru_3d,  gs, 2);
arv_3d  = create_field(arv_3d,  gs, 1);
cbc     = create_field(cbc,     gs, 3);
cor     = create_field(cor,     gs, 3);
d       = create_field(d,       gs, 3);
drx2d   = create_field(drx2d,   gs, 2);
dry2d   = create_field(dry2d,   gs, 1);
dt      = create_field(dt,      gs, 3);
fsm     = create_field(fsm,     gs, 3);
dum     = create_field(dum,     gs, 2);
dvm     = create_field(dvm,     gs, 1);

advx    = create_field(zeros(im,jm,kb),    gs, 2);
advy    = create_field(zeros(im,jm,kb),    gs, 1);

advua   = create_field(zeros(im,jm,1),   gs, 2);
advva   = create_field(zeros(im,jm,1),   gs, 1);
curv2d  = create_field(zeros(im,jm,1),   gs, 3);


drhox   = create_field(zeros(im,jm,kb),   gs, 2);
drhoy   = create_field(zeros(im,jm,kb),   gs, 1);

e_atmos = create_field(e_atmos, gs, 3);
egb     = create_field(egb,     gs, 3);
egf     = create_field(egf,     gs, 3);
el      = create_field(el,      gs, 3);
elb     = create_field(elb,     gs, 3);
elf     = create_field(elf,     gs, 3);
et      = create_field(et,      gs, 3);
etb     = create_field(etb,     gs, 3);
etf     = create_field(etf,     gs, 3);
fluxua  = create_field(fluxua,  gs, 2);
fluxva  = create_field(fluxva,  gs, 1);
h       = create_field(h,       gs, 3);
h_3d    = create_field(h_3d,    gs, 3);
psi     = create_field(psi,     gs, 3);
rot     = create_field(rot,     gs, 3);

swrad   = create_field(swrad,   gs, 7);
vfluxb  = create_field(vfluxb,  gs, 3);
vfluxf  = create_field(vfluxf,  gs, 3);

tsurf   = create_field(tsurf,   gs, 3);
ssurf   = create_field(ssurf,   gs, 3);
ua      = create_field(ua,      gs, 2);
uab     = create_field(uab,     gs, 2);
uaf     = create_field(uaf,     gs, 2);
utb     = create_field(utb,     gs, 2);
utf     = create_field(utf,     gs, 2);
va      = create_field(va,      gs, 1);
vab     = create_field(vab,     gs, 1);
vaf     = create_field(vaf,     gs, 1);
vtb     = create_field(vtb,     gs, 1);
vtf     = create_field(vtf,     gs, 1);

wssurf  = create_field(wssurf,  gs, 3);
wtsurf  = create_field(wtsurf,  gs, 3);
wubot   = create_field(wubot,   gs, 2);
wusurf  = create_field(wusurf,  gs, 2);
wvbot   = create_field(wvbot,   gs, 1);
wvsurf  = create_field(wvsurf,  gs, 1);
wubot1  = create_field(wubot1,  gs, 2);
wvbot1  = create_field(wubot1,  gs, 1);
aam     = create_field(aam,     gs, 3);

a       = create_field(a,       gs, 3);
c       = create_field(c,       gs, 3);

dtef    = create_field(dtef,    gs, 3);
ee      = create_field(ee,      gs, 3);
gg      = create_field(gg,      gs, 3);
kh      = create_field(kh,      gs, 7);
km      = create_field(km,      gs, 7);
kq      = create_field(kq,      gs, 7);
l       = create_field(l,       gs, 7);
q2b     = create_field(q2b,     gs, 7);
q2      = create_field(q2,      gs, 7);
q2lb    = create_field(q2lb,    gs, 7);
q2l     = create_field(q2l,     gs, 7);
rho     = create_field(rho,     gs, 3);
rmean   = create_field(rmean,   gs, 3);
sb      = create_field(sb,      gs, 3);
sclim   = create_field(sclim,   gs, 3);
s       = create_field(s,       gs, 3);
t       = create_field(t,       gs, 3);
tb      = create_field(tb,      gs, 3);
tclim   = create_field(tclim,   gs, 3);
ub      = create_field(ub,      gs, 2);
uf      = create_field(uf,      gs, 2);
u       = create_field(u,       gs, 2);
vb      = create_field(vb,      gs, 1);
vf      = create_field(vf,      gs, 1);
v       = create_field(v,       gs, 1);
w       = create_field(w,       gs, 7);
zflux   = create_field(zflux,   gs, 3);

ele     = create_field(ele,     gs, 3);
eln     = create_field(eln,     gs, 3);
els     = create_field(els,     gs, 3);
elw     = create_field(elw,     gs, 3);
sbe     = create_field(sbe,     gs, 3);
sbn     = create_field(sbn,     gs, 3);
sbs     = create_field(sbs,     gs, 3);
sbw     = create_field(sbw,     gs, 3);
tbe     = create_field(tbe,     gs, 3);
tbn     = create_field(tbn,     gs, 3);
tbs     = create_field(tbs,     gs, 3);
tbw     = create_field(tbw,     gs, 3);
uabe    = create_field(uabe,    gs, 3);
uabw    = create_field(uabw,    gs, 3);
ube     = create_field(ube,     gs, 3);
ubw     = create_field(ubw,     gs, 3);
vabn    = create_field(vabn,    gs, 3);
vabs    = create_field(vabs,    gs, 3);
vbn     = create_field(vbn,     gs, 3);
vbs     = create_field(vbs,     gs, 3);
d_3d    = create_field(d_3d,    gs, 3);
dt_3d   = create_field(dt_3d,   gs, 3);