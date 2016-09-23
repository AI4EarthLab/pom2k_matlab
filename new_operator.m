function [OP_AXF1_XY, OP_AXB1_XY, OP_AYF1_XY, OP_AYB1_XY, ...
          OP_AXF2_XY, OP_AXB2_XY, OP_AYF2_XY, OP_AYB2_XY, ...        
          OP_AXF1_XZ, OP_AXB1_XZ, OP_AZF1_XZ, OP_AZB1_XZ, ...
          OP_AXF2_XZ, OP_AXB2_XZ, OP_AZF2_XZ, OP_AZB2_XZ, ...
          OP_AYF1_YZ, OP_AYB1_YZ, OP_AZF1_YZ, OP_AZB1_YZ, ...
          OP_AYF2_YZ, OP_AYB2_YZ, OP_AZF2_YZ, OP_AZB2_YZ, ...
          OP_DXF1_XY, OP_DXB1_XY, OP_DYF1_XY, OP_DYB1_XY, ...
          OP_DXF2_XY, OP_DXB2_XY, OP_DYF2_XY, OP_DYB2_XY, ...
          OP_DXF1_XZ, OP_DXB1_XZ, OP_DZF1_XZ, OP_DZB1_XZ, ...
          OP_DXF2_XZ, OP_DXB2_XZ, OP_DZF2_XZ, OP_DZB2_XZ, ...
          OP_DYF1_YZ, OP_DYB1_YZ, OP_DZF1_YZ, OP_DZB1_YZ, ...
          OP_DYF2_YZ, OP_DYB2_YZ, OP_DZF2_YZ, OP_DZB2_YZ, ...
          OP_L_XY,   OP_L_XZ,   OP_L_YZ, ...
          OP_R_XY,   OP_R_XZ,   OP_R_YZ, ...
          OP_SUMX_XY, OP_SUMY_XY, ...
          OP_SUMX_XZ, OP_SUMZ_XZ, ...
          OP_SUMY_YZ, OP_SUMZ_YZ] = new_operator(m,n,k)
% Computing the basic averaging and differencing operators to smplify
% code.

% OP_AXF_XY: get the x average operator with forward method,    avg_f^x(i,j)=(f(i,j)+f(i+1,j))/2
% OP_AXB_XY: get the x average operator with backward method,   avg_f^x(i,j)=(f(i,j)+f(i-1,j))/2
% OP_AYF_XY: get the y average operator with forward method,    avg_f^y(i,j)=(f(i,j)+f(i,j+1))/2
% OP_AYB_XY: get the y average operator with backward method,   avg_f^y(i,j)=(f(i,j)+f(i,j-1))/2
% OP_DXF_XY: get the x difference operator with forward method, delta_x f(i,j)=(f(i+1,j)-f(i,  j))/2
% OP_DXB_XY: get the x difference operator with backward method,delta_x f(i,j)=(f(i,j  )-f(i-1,j))/2
% OP_DYF_XY: get the y difference operator with forward method, delta_y f(i,j)=(f(i,j+1)-f(i,  j))/2
% OP_DYB_XY: get the y difference operator with backward method,delta_y f(i,j)=(f(i,j  )-f(i,j-1))/2
%
% OP_AXF_XZ: get the x average operator with forward method,    avg_f^x(i,k)=(f(i,k)+f(i+1,j))/2
% OP_AXB_XZ: get the x average operator with backward method,   avg_f^x(i,k)=(f(i,k)+f(i-1,j))/2
% OP_AYF_XZ: get the y average operator with forward method,    avg_f^y(i,k)=(f(i,k)+f(i,j+1))/2
% OP_AYB_XZ: get the y average operator with backward method,   avg_f^y(i,k)=(f(i,k)+f(i,j-1))/2
% OP_DXF_XZ: get the x difference operator with forward method, delta_x f(i,k)=(f(i+1,k)-f(i,  k))/2
% OP_DXB_XZ: get the x difference operator with backward method,delta_x f(i,k)=(f(i,  k)-f(i-1,k))/2
% OP_DYF_XZ: get the y difference operator with forward method, delta_y f(i,k)=(f(i,k+1)-f(i,  k))/2
% OP_DYB_XZ: get the y difference operator with backward method,delta_y f(i,k)=(f(i,k  )-f(i,k-1))/2
%
% OP_AXF_XY: get the x average operator with forward method,    avg_f^x(j,k)=(f(j,k)+f(j+1,k))/2
% OP_AXB_XY: get the x average operator with backward method,   avg_f^x(j,k)=(f(j,k)+f(j-1,k))/2
% OP_AYF_XY: get the y average operator with forward method,    avg_f^y(j,k)=(f(j,k)+f(j,k+1))/2
% OP_AYB_XY: get the y average operator with backward method,   avg_f^y(j,k)=(f(j,k)+f(j,k-1))/2
% OP_DXF_XY: get the x difference operator with forward method, delta_x f(j,k)=(f(j+1,k)-f(j,  k))/2
% OP_DXB_XY: get the x difference operator with backward method,delta_x f(j,k)=(f(j,  k)-f(j-1,k))/2
% OP_DYF_XY: get the y difference operator with forward method, delta_y f(j,k)=(f(j,k+1)-f(j,  k))/2
% OP_DYB_XY: get the y difference operator with backward method,delta_y f(j,k)=(f(j,  k)-f(j,k-1))/2
%


% Suppose:
% X= [X11 X12 X13 X14 X15]      Y= [Y11 Y12 Y13 Y14 Y15]    Z= [Z11 Z12 Z13 Z14 Z15 Z16]
%    [X21 X22 X23 X24 X25]         [Y21 Y22 Y23 Y24 Y25]       [Z21 Z22 Z23 Z24 Z25 Z26]
%    [X31 X32 X33 X34 X35]         [Y31 Y32 Y33 Y34 Y35]       [Z31 Z32 Z33 Z34 Z35 Z36]
%    [X41 X42 X43 X44 X45]         [Y41 Y42 Y43 Y44 Y45]       [Z41 Z42 Z43 Z44 Z45 Z46]
%    [X51 X52 X53 X54 X55]         [Y51 Y52 Y53 Y54 Y55]       [Z51 Z52 Z53 Z54 Z55 Z56]
%    [X61 X62 X63 X64 X65]         [Y61 Y62 Y63 Y64 Y65]       [Z61 Z62 Z63 Z64 Z65 Z66]
%    [X71 X72 X73 X74 X75]         [Y71 Y72 Y73 Y74 Y75]       [Z71 Z72 Z73 Z74 Z75 Z76]


%%%%%%%%%%%%%%%%%%%%%%%% A operator in XY plate %%%%%%%%%%%%%%%%%%
% OP_AXF_XY1= 0.5*[1  1  0  0  0  0  0]  OP_AXF_XY1*X= 0.5*[ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%                 [0  1  1  0  0  0  0]                    [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%                 [0  0  1  1  0  0  0]                    [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%                 [0  0  0  1  1  0  0]                    [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%                 [0  0  0  0  1  1  0]                    [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%                 [0  0  0  0  0  1  1]                    [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]
%                 [0  0  0  0  0  0  0]                    [       0        0        0        0        0 ]
L1= zeros(m,m); L1(1:m-1, 1:m-1) = eye(m-1);
L2= zeros(m,m); L2(1:m-1, 2:m  ) = eye(m-1);
OP_AXF1_XY=0.5*(L1+L2);

% OP_AXB1_XY= 0.5*[0  0  0  0  0  0  0]  OP_AXB1_XY*X= 0.5*[       0        0        0        0        0 ]
%                 [1  1  0  0  0  0  0]                    [ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%                 [0  1  1  0  0  0  0]                    [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%                 [0  0  1  1  0  0  0]                    [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%                 [0  0  0  1  1  0  0]                    [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%                 [0  0  0  0  1  1  0]                    [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%                 [0  0  0  0  0  1  1]                    [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]
L1= zeros(m,m); L1(2:m, 1:m-1) = eye(m-1);
L2= zeros(m,m); L2(2:m, 2:m  ) = eye(m-1);
OP_AXB1_XY=0.5*(L1+L2);

% OP_AYF1_XY= 0.5*[1  0  0  0  0]        Y*OP_AYF1_XY= 0.5*[ Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15  0 ]
%                 [1  1  0  0  0]                          [ Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25  0 ]
%                 [0  1  1  0  0]                          [ Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35  0 ]
%                 [0  0  1  1  0]                          [ Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45  0 ]
%                 [0  0  0  1  0]                          [ Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55  0 ]
%                                                          [ Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65  0 ]
%                                                          [ Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75  0 ]
R1= zeros(n,n); R1(1:n-1, 1:n-1) = eye(n-1);
R2= zeros(n,n); R2(2:n  , 1:n-1) = eye(n-1);
OP_AYF1_XY = 0.5*(R1+R2);

% OP_AYB1_XY= 0.5*[0  1  0  0  0]        Y*OP_AYB1_XY= 0.5*[ 0  Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15 ]
%                 [0  1  1  0  0]                          [ 0  Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25 ]
%                 [0  0  1  1  0]                          [ 0  Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35 ]
%                 [0  0  0  1  1]                          [ 0  Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45 ]
%                 [0  0  0  0  1]                          [ 0  Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55 ]
%                                                          [ 0  Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65 ]
%                                                          [ 0  Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75 ]
R1= zeros(n,n); R1(1:n-1, 2:n) = eye(n-1);
R2= zeros(n,n); R2(2:n  , 2:n) = eye(n-1);
OP_AYB1_XY = 0.5*(R1+R2);

% OP_AXF2_XY= 0.5*[0  0  0  0  0  0  0]  OP_AXF2_XY*X= 0.5*[       0        0        0        0        0 ]
%                 [0  1  1  0  0  0  0]                    [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%                 [0  0  1  1  0  0  0]                    [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%                 [0  0  0  1  1  0  0]                    [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%                 [0  0  0  0  1  1  0]                    [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%                 [0  0  0  0  0  1  1]                    [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]
%                 [0  0  0  0  0  0  0]                    [       0        0        0        0        0 ]
L1= zeros(m,m); L1(2:m-1, 2:m-1) = eye(m-2);
L2= zeros(m,m); L2(2:m-1, 3:m  ) = eye(m-2);
OP_AXF2_XY=0.5*(L1+L2);

% OP_AXB2_XY= 0.5*[0  0  0  0  0  0  0]  OP_AXB2_XY*X= 0.5*[       0        0        0        0        0 ]
%                 [1  1  0  0  0  0  0]                    [ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%                 [0  1  1  0  0  0  0]                    [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%                 [0  0  1  1  0  0  0]                    [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%                 [0  0  0  1  1  0  0]                    [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%                 [0  0  0  0  1  1  0]                    [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%                 [0  0  0  0  0  0  0]                    [       0        0        0        0        0 ]
L1= zeros(m,m); L1(2:m-1, 1:m-2) = eye(m-2);
L2= zeros(m,m); L2(2:m-1, 2:m-1) = eye(m-2);
OP_AXB2_XY=0.5*(L1+L2);

% OP_AYF2_XY= 0.5*[0  0  0  0  0]        Y*OP_AYF2_XY= 0.5*[ 0  Y12+Y13  Y13+Y14  Y14+Y15  0 ]
%                 [0  1  0  0  0]                          [ 0  Y22+Y23  Y23+Y24  Y24+Y25  0 ]
%                 [0  1  1  0  0]                          [ 0  Y32+Y33  Y33+Y34  Y34+Y35  0 ]
%                 [0  0  1  1  0]                          [ 0  Y42+Y43  Y43+Y44  Y44+Y45  0 ]
%                 [0  0  0  1  0]                          [ 0  Y52+Y53  Y53+Y54  Y54+Y55  0 ]
%                                                          [ 0  Y62+Y63  Y63+Y64  Y64+Y65  0 ]
%                                                          [ 0  Y72+Y73  Y73+Y74  Y74+Y75  0 ]
R1= zeros(n,n); R1(2:n-1, 2:n-1) = eye(n-2);
R2= zeros(n,n); R2(3:n  , 2:n-1) = eye(n-2);
OP_AYF2_XY = 0.5*(R1+R2);

% OP_AYB2_XY= 0.5*[0  1  0  0  0]        Y*OP_AYB2_XY= 0.5*[ 0  Y11+Y12  Y12+Y13  Y13+Y14  0 ]
%                 [0  1  1  0  0]                          [ 0  Y21+Y22  Y22+Y23  Y23+Y24  0 ]
%                 [0  0  1  1  0]                          [ 0  Y31+Y32  Y32+Y33  Y33+Y34  0 ]
%                 [0  0  0  1  0]                          [ 0  Y41+Y42  Y42+Y43  Y43+Y44  0 ]
%                 [0  0  0  0  0]                          [ 0  Y51+Y52  Y52+Y53  Y53+Y54  0 ]
%                                                          [ 0  Y61+Y62  Y62+Y63  Y63+Y64  0 ]
%                                                          [ 0  Y71+Y72  Y72+Y73  Y73+Y74  0 ]
R1= zeros(n,n); R1(1:n-2, 2:n-1) = eye(n-2);
R2= zeros(n,n); R2(2:n-1, 2:n-1) = eye(n-2);
OP_AYB2_XY = 0.5*(R1+R2);


%%%%%%%%%%%%%%%%%%%%%%%% A operator in XZ plate %%%%%%%%%%%%%%%%%%
L1= zeros(m,m); L1(1:m-1, 1:m-1) = eye(m-1);
L2= zeros(m,m); L2(1:m-1, 2:m  ) = eye(m-1);
OP_AXF1_XZ=0.5*(L1+L2);

L1= zeros(m,m); L1(2:m, 1:m-1) = eye(m-1);
L2= zeros(m,m); L2(2:m, 2:m  ) = eye(m-1);
OP_AXB1_XZ=0.5*(L1+L2);

R1= zeros(k,k); R1(1:k-1, 1:k-1) = eye(k-1);
R2= zeros(k,k); R2(2:k  , 1:k-1) = eye(k-1);
OP_AZF1_XZ = 0.5*(R1+R2);

R1= zeros(k,k); R1(1:k-1, 2:k) = eye(k-1);
R2= zeros(k,k); R2(2:k  , 2:k) = eye(k-1);
OP_AZB1_XZ = 0.5*(R1+R2);

L1= zeros(m,m); L1(2:m-1, 2:m-1) = eye(m-2);
L2= zeros(m,m); L2(2:m-1, 3:m  ) = eye(m-2);
OP_AXF2_XZ=0.5*(L1+L2);

L1= zeros(m,m); L1(2:m-1, 1:m-2) = eye(m-2);
L2= zeros(m,m); L2(2:m-1, 2:m-1) = eye(m-2);
OP_AXB2_XZ=0.5*(L1+L2);

R1= zeros(k,k); R1(2:k-1, 2:k-1) = eye(k-2);
R2= zeros(k,k); R2(3:k  , 2:k-1) = eye(k-2);
OP_AZF2_XZ = 0.5*(R1+R2);

R1= zeros(k,k); R1(1:k-2, 2:k-1) = eye(k-2);
R2= zeros(k,k); R2(2:k-1, 2:k-1) = eye(k-2);
OP_AZB2_XZ = 0.5*(R1+R2);

%%%%%%%%%%%%%%%%%%%%%%%% A operator in YZ plate %%%%%%%%%%%%%%%%%%
L1= zeros(n,n); L1(1:n-1, 1:n-1) = eye(n-1);
L2= zeros(n,n); L2(1:n-1, 2:n  ) = eye(n-1);
OP_AYF1_YZ=0.5*(L1+L2);

L1= zeros(n,n); L1(2:n, 1:n-1) = eye(n-1);
L2= zeros(n,n); L2(2:n, 2:n  ) = eye(n-1);
OP_AYB1_YZ=0.5*(L1+L2);

R1= zeros(k,k); R1(1:k-1, 1:k-1) = eye(k-1);
R2= zeros(k,k); R2(2:k  , 1:k-1) = eye(k-1);
OP_AZF1_YZ = 0.5*(R1+R2);

R1= zeros(k,k); R1(1:k-1, 2:k) = eye(k-1);
R2= zeros(k,k); R2(2:k  , 2:k) = eye(k-1);
OP_AZB1_YZ = 0.5*(R1+R2);

L1= zeros(n,n); L1(2:n-1, 2:n-1) = eye(n-2);
L2= zeros(n,n); L2(2:n-1, 3:n  ) = eye(n-2);
OP_AYF2_YZ=0.5*(L1+L2);


L1= zeros(n,n); L1(2:n-1, 1:n-2) = eye(n-2);
L2= zeros(n,n); L2(2:n-1, 2:n-1) = eye(n-2);
OP_AYB2_YZ=0.5*(L1+L2);

R1= zeros(k,k); R1(2:k-1, 2:k-1) = eye(k-2);
R2= zeros(k,k); R2(3:k  , 2:k-1) = eye(k-2);
OP_AZF2_YZ = 0.5*(R1+R2);


R1= zeros(k,k); R1(1:k-2, 2:k-1) = eye(k-2);
R2= zeros(k,k); R2(2:k-1 ,2:k-1) = eye(k-2);
OP_AZB2_YZ = 0.5*(R1+R2);


%%%%%%%%%%%%%%%%%%%%%%%% D operator in XY plate %%%%%%%%%%%%%%%%%%
% OP_DXF1_XY=[ -1  1  0  0  0  0  0]    OP_DXF1_XY*X =[ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%            [  0 -1  1  0  0  0  0]                  [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%            [  0  0 -1  1  0  0  0]                  [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%            [  0  0  0 -1  1  0  0]                  [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%            [  0  0  0  0 -1  1  0]                  [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%            [  0  0  0  0  0 -1  1]                  [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
%            [  0  0  0  0  0  0  0]                  [       0        0         0         0         0 ]
L1= zeros(m,m); L1(1:m-1, 1:m-1) = -eye(m-1);
L2= zeros(m,m); L2(1:m-1, 2:m  ) =  eye(m-1);
OP_DXF1_XY= L1+L2;


% OP_DXB1_XY=[  0  0  0  0  0  0  0]    OP_DXB1_XY*X =[       0        0         0         0         0 ]
%            [ -1  1  0  0  0  0  0]                  [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%            [  0 -1  1  0  0  0  0]                  [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%            [  0  0 -1  1  0  0  0]                  [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%            [  0  0  0 -1  1  0  0]                  [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%            [  0  0  0  0 -1  1  0]                  [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%            [  0  0  0  0  0 -1  1]                  [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
L1= zeros(m,m); L1(2:m, 1:m-1) = -eye(m-1);
L2= zeros(m,m); L2(2:m, 2:m  ) =  eye(m-1);
OP_DXB1_XY= L1+L2;

% OP_DYF1_XY=[ -1  0  0  0  0 ]         Y*OP_DYF1_XY =[ Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14  0 ]
%            [  1 -1  0  0  0 ]                       [ Y22-Y21   Y23-Y22   Y24-Y23   Y25-Y24  0 ]
%            [  0  1 -1  0  0 ]                       [ Y22-Y21   Y33-Y32   Y34-Y33   Y35-Y34  0 ]
%            [  0  0  1 -1  0 ]                       [ Y22-Y21   Y43-Y42   Y44-Y43   Y45-Y44  0 ]
%            [  0  0  0  1  0 ]                       [ Y22-Y21   Y53-Y52   Y54-Y53   Y55-Y54  0 ]
%                                                     [ Y22-Y21   Y63-Y62   Y64-Y63   Y65-Y64  0 ]
%                                                     [ Y22-Y21   Y73-Y72   Y74-Y73   Y75-Y74  0 ]
R1= zeros(n,n); R1(1:n-1, 1:n-1) = -eye(n-1);
R2= zeros(n,n); R2(2:n  , 1:n-1) =  eye(n-1);
OP_DYF1_XY= R1+R2;

% OP_DYB1_XY=[  0 -1  0  0  0 ]         Y*OP_DYB1_XY =[ 0  Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14 ]
%            [  0  1 -1  0  0 ]                       [ 0  Y22-Y21   Y23-Y22   Y24-Y23   Y15-Y14 ]
%            [  0  0  1 -1  0 ]                       [ 0  Y32-Y31   Y33-Y32   Y34-Y33   Y35-Y34 ]
%            [  0  0  0  1 -1 ]                       [ 0  Y42-Y41   Y43-Y42   Y44-Y43   Y45-Y44 ]
%            [  0  0  0  0  1 ]                       [ 0  Y52-Y51   Y53-Y52   Y54-Y53   Y55-Y54 ]
%                                                     [ 0  Y62-Y61   Y63-Y62   Y64-Y63   Y65-Y64 ]
%                                                     [ 0  Y72-Y71   Y73-Y72   Y74-Y73   Y75-Y74 ]
R1= zeros(n,n); R1(1:n-1, 2:n) = -eye(n-1);
R2= zeros(n,n); R2(2:n  , 2:n) =  eye(n-1);
OP_DYB1_XY= R1+R2;



% OP_DXF2_XY=[  0  0  0  0  0  0  0]    OP_DXF2_XY*X =[       0        0         0         0         0 ]
%            [  0 -1  1  0  0  0  0]                  [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%            [  0  0 -1  1  0  0  0]                  [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%            [  0  0  0 -1  1  0  0]                  [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%            [  0  0  0  0 -1  1  0]                  [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%            [  0  0  0  0  0 -1  1]                  [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
%            [  0  0  0  0  0  0  0]                  [       0        0         0         0         0 ]
L1= zeros(m,m); L1(2:m-1, 2:m-1) = -eye(m-2);
L2= zeros(m,m); L2(2:m-1, 3:m  ) =  eye(m-2);
OP_DXF2_XY= L1+L2;


% OP_DXB2_XY=[  0  0  0  0  0  0  0]    OP_DXB2_XY*X =[       0        0         0         0         0 ]
%            [ -1  1  0  0  0  0  0]                  [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%            [  0 -1  1  0  0  0  0]                  [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%            [  0  0 -1  1  0  0  0]                  [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%            [  0  0  0 -1  1  0  0]                  [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%            [  0  0  0  0 -1  1  0]                  [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%            [  0  0  0  0  0  0  0]                  [       0        0         0         0         0 ]
L1= zeros(m,m); L1(2:m-1, 1:m-2) = -eye(m-2);
L2= zeros(m,m); L2(2:m-1, 2:m-1) =  eye(m-2);
OP_DXB2_XY= L1+L2;

% OP_DYF2_XY=[  0  0  0  0  0 ]         Y*OP_DYF2_XY =[ 0   Y13-Y12   Y14-Y13   Y15-Y14  0 ]
%            [  0 -1  0  0  0 ]                       [ 0   Y23-Y22   Y24-Y23   Y25-Y24  0 ]
%            [  0  1 -1  0  0 ]                       [ 0   Y33-Y32   Y34-Y33   Y35-Y34  0 ]
%            [  0  0  1 -1  0 ]                       [ 0   Y43-Y42   Y44-Y43   Y45-Y44  0 ]
%            [  0  0  0  1  0 ]                       [ 0   Y53-Y52   Y54-Y53   Y55-Y54  0 ]
%                                                     [ 0   Y63-Y62   Y64-Y63   Y65-Y64  0 ]
%                                                     [ 0   Y73-Y72   Y74-Y73   Y75-Y74  0 ]
R1= zeros(n,n); R1(2:n-1, 2:n-1) = -eye(n-2);
R2= zeros(n,n); R2(3:n  , 2:n-1) =  eye(n-2);
OP_DYF2_XY= R1+R2;

% OP_DYB2_XY=[  0 -1  0  0  0 ]         Y*OP_DYB2_XY =[ 0  Y12-Y11   Y13-Y12   Y14-Y13   0 ]
%            [  0  1 -1  0  0 ]                       [ 0  Y22-Y21   Y23-Y22   Y24-Y23   0 ]
%            [  0  0  1 -1  0 ]                       [ 0  Y32-Y31   Y33-Y32   Y34-Y33   0 ]
%            [  0  0  0  1  0 ]                       [ 0  Y42-Y41   Y43-Y42   Y44-Y43   0 ]
%            [  0  0  0  0  0 ]                       [ 0  Y52-Y51   Y53-Y52   Y54-Y53   0 ]
%                                                     [ 0  Y62-Y61   Y63-Y62   Y64-Y63   0 ]
%                                                     [ 0  Y72-Y71   Y73-Y72   Y74-Y73   0 ]
R1= zeros(n,n); R1(1:n-2, 2:n-1) = -eye(n-2);
R2= zeros(n,n); R2(2:n-1, 2:n-1) =  eye(n-2);
OP_DYB2_XY= R1+R2;

%%%%%%%%%%%%%%%%%%%%%%%% D operator in XZ plate %%%%%%%%%%%%%%%%%%
L1= zeros(m,m); L1(1:m-1, 1:m-1) = -eye(m-1);
L2= zeros(m,m); L2(1:m-1, 2:m  ) =  eye(m-1);
OP_DXF1_XZ= L1+L2;

L1= zeros(m,m); L1(2:m, 1:m-1) = -eye(m-1);
L2= zeros(m,m); L2(2:m, 2:m  ) =  eye(m-1);
OP_DXB1_XZ= L1+L2;

R1= zeros(k,k); R1(1:k-1, 1:k-1) = -eye(k-1);
R2= zeros(k,k); R2(2:k  , 1:k-1) =  eye(k-1);
OP_DZF1_XZ= R1+R2;

R1= zeros(k,k); R1(1:k-1, 2:k ) = -eye(k-1);
R2= zeros(k,k); R2(2:k  , 2:k ) =  eye(k-1);
OP_DZB1_XZ= R1+R2;

L1= zeros(m,m); L1(2:m-1, 2:m-1) = -eye(m-2);
L2= zeros(m,m); L2(2:m-1, 3:m  ) =  eye(m-2);
OP_DXF2_XZ= L1+L2;

L1= zeros(m,m); L1(2:m-1, 1:m-2) = -eye(m-2);
L2= zeros(m,m); L2(2:m-1, 2:m-1) =  eye(m-2);
OP_DXB2_XZ= L1+L2;

R1= zeros(k,k); R1(2:k-1, 2:k-1) = -eye(k-2);
R2= zeros(k,k); R2(3:k  , 2:k-1) =  eye(k-2);
OP_DZF2_XZ= R1+R2;

R1= zeros(k,k); R1(1:k-2, 2:k-1) = -eye(k-2);
R2= zeros(k,k); R2(2:k-1, 2:k-1) =  eye(k-2);
OP_DZB2_XZ= R1+R2;



%%%%%%%%%%%%%%%%%%%%%%%% D operator in YZ plate %%%%%%%%%%%%%%%%%%
L1= zeros(n,n); L1(1:n-1, 1:n-1) = -eye(n-1);
L2= zeros(n,n); L2(1:n-1, 2:n  ) =  eye(n-1);
OP_DYF1_YZ= L1+L2;

L1= zeros(n,n); L1(2:n, 1:n-1) = -eye(n-1);
L2= zeros(n,n); L2(2:n, 2:n  ) =  eye(n-1);
OP_DYB1_YZ= L1+L2;

R1= zeros(k,k); R1(1:k-1, 1:k-1) = -eye(k-1);
R2= zeros(k,k); R2(2:k  , 1:k-1) =  eye(k-1);
OP_DZF1_YZ= R1+R2;

R1= zeros(k,k); R1(1:k-1, 2:k) = -eye(k-1);
R2= zeros(k,k); R2(2:k  , 2:k) =  eye(k-1);
OP_DZB1_YZ= R1+R2;

L1= zeros(n,n); L1(2:n-1, 2:n-1) = -eye(n-2);
L2= zeros(n,n); L2(2:n-1, 3:n  ) =  eye(n-2);
OP_DYF2_YZ= L1+L2;

L1= zeros(n,n); L1(2:n-1, 1:n-2) = -eye(n-2);
L2= zeros(n,n); L2(2:n-1, 2:n-1) =  eye(n-2);
OP_DYB2_YZ= L1+L2;

R1= zeros(k,k); R1(2:k-1, 2:k-1) = -eye(k-2);
R2= zeros(k,k); R2(3:k  , 2:k-1) =  eye(k-2);
OP_DZF2_YZ= R1+R2;

R1= zeros(k,k); R1(1:k-2, 2:k-1) = -eye(k-2);
R2= zeros(k,k); R2(2:k-1, 2:k-1) =  eye(k-2);
OP_DZB2_YZ= R1+R2;

%%%%%%%%%%%%%%%%%%%%%%%% OP_L and OP_R operator %%%%%%%%%%%%%%%%%%
L=eye(m,m); L(1,1)=0;L(m,m)=0;
R=eye(n,n); R(1,1)=0;R(n,n)=0;
OP_L_XY = L;
OP_R_XY = R;

L=eye(m,m); L(1,1)=0;L(m,m)=0;
R=eye(k,k); R(1,1)=0;R(k,k)=0;
OP_L_XZ = L;
OP_R_XZ = R;

L=eye(n,n); L(1,1)=0;L(n,n)=0;
R=eye(k,k); R(1,1)=0;R(k,k)=0;
OP_L_YZ = L;
OP_R_YZ = R;

%%%%%%%%%%%%%%%%%%%%%%%% OP_SUM operator %%%%%%%%%%%%%%%%%%
% Define sum operator along with X and Y directions. 
%   L=[1 0 0 0 0 0 0]      R=[1 1 1 1 0]
%     [1 1 0 0 0 0 0]        [0 1 1 1 0] 
%     [1 1 1 0 0 0 0]        [0 0 1 1 0]
%     [1 1 1 1 0 0 0]        [0 0 0 1 0]
%     [1 1 1 1 1 0 0]        [0 0 0 0 0]
%     [1 1 1 1 1 1 0]        
%     [0 0 0 0 0 0 0]        
L=tril(ones(m,m)); L(m,:)=0;
R=triu(ones(n,n)); R(:,n)=0;
OP_SUMX_XY= L;
OP_SUMY_XY= R;

L=tril(ones(m,m)); L(m,:)=0;
R=triu(ones(k,k)); R(:,k)=0;
OP_SUMX_XZ= L;
OP_SUMZ_XZ= R;


L=tril(ones(n,n)); L(n,:)=0;
R=triu(ones(k,k)); R(:,k)=0;
OP_SUMY_YZ= L;
OP_SUMZ_YZ= R;

end

function test_a()
F=OP_AXB_XY*X;
end



