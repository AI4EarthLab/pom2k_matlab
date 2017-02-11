m = 5;
n = 5;
k = 3;

%***********************************
dx=0.11*ones(m,n);  %dx_2d
dy=0.12*ones(m,n);  %dy_2d
dz=0.13*ones(1,k);  %dz_2d

%***********************************
DATA_U   = rand(m, n, k);
DATA_UB  = rand(m, n, k);
DATA_D   = rand(m, n, k);
DATA_V   = rand(m, n, k);
DATA_VB  = rand(m, n, k);
DATA_AAM = rand(m, n, k);

%gs = create_grid(m, n, k, 'C');

%matlab does not support modifying an input parameter inside call
gs = init_grid('C', dx, dy, dz);

D   = create_field(DATA_D,   gs, 3);
AAM = create_field(DATA_AAM, gs, 3);
V   = create_field(DATA_V,   gs, 1);
VB  = create_field(DATA_VB,  gs, 1);
U   = create_field(DATA_U,   gs, 2);
UB  = create_field(DATA_UB,  gs, 2);


%examples
xflux1 = DXB(AXF(AXB(D) .* U) .* AXF(U))  - DXB(D .* AAM .* 2.0 .* DXF(UB));
yflux1 = DYB(AYF(AYB(D) .* V) .* AYF(V))  - DYB(D .* AAM .* 2.0 .* DYF(VB));

