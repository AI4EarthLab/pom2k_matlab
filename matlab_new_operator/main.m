m = 5;
n = 5;
k = 3;

%***********************************
v1  = repmat(0.11, m, 1);
v11 = ones(m, 1);

v2 = repmat(0.12, n, 1);
v21 = ones(n, 1);

v3 = repmat(0.12, k, 1);
v31 = ones(k, 1);

dx = meshgrid(v1, v21, v31); % dx_3d 
dy = meshgrid(v1, v21, v31); % dy_3d
dz = meshgrid(v1, v21, v31); % dz_3d

%***********************************
DATA_U   = rand(m, n, k);
DATA_UB  = rand(m, n, k);
DATA_D   = rand(m, n, k);
DATA_V   = rand(m, n, k);
DATA_VB  = rand(m, n, k);
DATA_AAM = rand(m, n, k);

gs = create_grid(m, n, k, 'C');

%matlab does not support modifying an input parameter inside call
gs = init_grid(gs, dx, dy, dz);

D   = create_field(DATA_D,   gs, 3);
AAM = create_field(DATA_AAM, gs, 3);
V   = create_field(DATA_V,   gs, 1);
VB  = create_field(DATA_VB,  gs, 1);
U   = create_field(DATA_U,   gs, 2);
UB  = create_field(DATA_UB,  gs, 2);


%examples
xflux1 = DXB(AXF(AXB(D) .* U) .* AXF(U))  - DXB(D .* AAM .* 2.0 .* DXF(UB));
yflux1 = DYB(AYF(AYB(D) .* V) .* AYF(V))  - DYB(D .* AAM .* 2.0 .* DYF(VB));

