m = 5;
n = 5;
k = 3;

%***********************************
v1  = repmat(2, m, 1);
v11 = ones(m, 1);

v2 = repmat(3, n, 1);
v21 = ones(n, 1);

v3 = repmat(4, k, 1);
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
gs = init_grid(gs, dx, dy, dz);

D   = create_field(DATA_D,   gs, 3);
AAM = create_field(DATA_AAM, gs, 3);
V   = create_field(DATA_V,   gs, 1);
VB  = create_field(DATA_VB,  gs, 1);
U   = create_field(DATA_U,   gs, 2);
UB  = create_field(DATA_UB,  gs, 2);

T   = create_field(DATA_D,   gs, 7);

%examples
xflux1 = DXB(AXF(AXB(D) .* U) .* AXF(U))  - DXB(D .* AAM .* 2.0 .* DXF(UB));
yflux1 = DYB(AYF(AYB(D) .* V) .* AYF(V))  - DYB(D .* AAM .* 2.0 .* DYF(VB));

xflux3 = DXB(T);

X0=create_field(DATA_D,   gs, 0); 
X1=create_field(DATA_D,   gs, 1); 
X2=create_field(DATA_D,   gs, 2); 
X3=create_field(DATA_D,   gs, 3); 
X4=create_field(DATA_D,   gs, 4); 
X5=create_field(DATA_D,   gs, 5); 
X6=create_field(DATA_D,   gs, 6); 
X7=create_field(DATA_D,   gs, 7); 

DXB(X0); DXB(X1); DXB(X2); DXB(X3); DXB(X4); DXB(X5); DXB(X6);
DXF(X0); DXF(X1); DXF(X2); DXF(X3); DXF(X4); DXF(X5); DXF(X6);

DYB(X0); DYB(X1); DYB(X2); DYB(X3); DYB(X4); DYB(X5); DYB(X6);
DYF(X0); DYF(X1); DYF(X2); DYF(X3); DYF(X4); DYF(X5); DYF(X6);


