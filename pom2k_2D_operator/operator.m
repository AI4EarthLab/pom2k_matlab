function [AXF, AXB, AYF, AYB, DXF, DXB, DYF, DYB]=operator(m,n)
% Computing the basic averaging and differencing operators to simplify
% code.

% AXF: get the x average operator with forward method,  avg_f^x(i,j)=(f(i,j)+f(i+1,j))/2
% AXB: get the x average operator with backward method, avg_f^x(i,j)=(f(i,j)+f(i-1,j))/2
% AYF: get the y average operator with forward method, avg_f^y(i,j)=(f(i,j)+f(i,j+1))/2
% AYB: get the y average operator with backward method, avg_f^y(i,j)=(f(i,j)+f(i,j-1))/2

% DXF: get the x difference operator with forward method,  delta_x f(i,j)=(f(i+1,j)-f(i,j))/2
% DXB: get the x difference operator with backward method, delta_x f(i,j)=(f(i,j)-f(i-1,j))/2
% DYF: get the y difference operator with forward method,  delta_y f(i,j)=(f(i,j+1)-f(i,j))/2
% DYB: get the y difference operator with backward method, delta_y f(i,j)=(f(i,j)-f(i,j-1))/2

% X= [X11 X12 X13 X14 X15]      Y= [Y11 Y12 Y13 Y14 Y15]
%    [X21 X22 X23 X24 X25]         [Y21 Y22 Y23 Y24 Y25]
%    [X31 X32 X33 X34 X35]         [Y31 Y32 Y33 Y34 Y35]
%    [X41 X42 X43 X44 X45]         [Y41 Y42 Y43 Y44 Y45]
%    [X51 X52 X53 X54 X55]         [Y51 Y52 Y53 Y54 Y55]
%    [X61 X62 X63 X64 X65]         [Y61 Y62 Y63 Y64 Y65]
%    [X71 X72 X73 X74 X75]         [Y71 Y72 Y73 Y74 Y75]
%
% AXF= 0.5*[1  1  0  0  0  0  0]   AXF*X= 0.5*[ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%          [0  1  1  0  0  0  0]              [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%          [0  0  1  1  0  0  0]              [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%          [0  0  0  1  1  0  0]              [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%          [0  0  0  0  1  1  0]              [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%          [0  0  0  0  0  1  1]              [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]
%          [0  0  0  0  0  0  1]              [     X71      X72      X73      X74      X75 ]

L1= eye(m);
L2= zeros(m,m); L2(1:m-1,2:m) = eye(m-1);
%L3= zeros(m,m); L3(m,m)=-1;
AXF=0.5*(L1+L2);

% AXB= 0.5*[1  0  0  0  0  0  0]   AXB*X= 0.5*[     X11      X12      X13      X14      X15 ]
%          [1  1  0  0  0  0  0]              [ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%          [0  1  1  0  0  0  0]              [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%          [0  0  1  1  0  0  0]              [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%          [0  0  0  1  1  0  0]              [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%          [0  0  0  0  1  1  0]              [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%          [0  0  0  0  0  1  1]              [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]

L1= eye(m);
L2= zeros(m,m); L2(2:m, 1:m-1) = eye(m-1);
%L3= zeros(m,m); L3(1,1)=-1;
AXB=0.5*(L1+L2);

% AYF=0.5*[1  0  0  0  0]   Y*AYF*R=0.5*[ Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15  Y15 ]
%         [1  1  0  0  0]             [ Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25  Y25 ]
%         [0  1  1  0  0]             [ Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35  Y35 ]
%         [0  0  1  1  0]             [ Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45  Y45 ]
%         [0  0  0  1  1]             [ Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55  Y55 ]
%                                     [ Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65  Y65 ]
%                                     [ Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75  Y75 ]

R1= eye(n);
R2= zeros(n,n); R2(2:n,1:n-1)= eye(n-1);
%R3= zeros(n,n); R3(n,n)=-1;
AYF = 0.5*(R1+R2);

% AYB=0.5*[1  1  0  0  0]   Y*AYB=0.5*[ Y11  Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15 ]
%         [0  1  1  0  0]             [ Y12  Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25 ]
%         [0  0  1  1  0]             [ Y13  Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35 ]
%         [0  0  0  1  1]             [ Y14  Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45 ]
%         [0  0  0  0  1]             [ Y15  Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55 ]
%                                     [ Y16  Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65 ]
%                                     [ Y17  Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75 ]

R1= eye(n);
R2= zeros(n,n); R2(1:n-1,2:n)= eye(n-1);
%R3= zeros(n,n); R3(1,1)=-1;
AYB = 0.5*(R1+R2);

% DXF=[ -1  1  0  0  0  0  0]   DXF*X=[ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%     [  0 -1  1  0  0  0  0]         [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%     [  0  0 -1  1  0  0  0]         [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%     [  0  0  0 -1  1  0  0]         [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%     [  0  0  0  0 -1  1  0]         [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%     [  0  0  0  0  0 -1  1]         [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
%     [  0  0  0  0  0  0 -1]         [    -X71     -X72      -X73      -X74      -X75 ]
									  
L1= -eye(m);
L2= zeros(m,m); L2(1:m-1,2:m) = eye(m-1);
%L3= zeros(m,m); L3(m,m)=1;
DXF=(L1+L2);

% DXB=[  1  0  0  0  0  0  0]   DXB*X=[     X11      X12       X13       X14       X15 ]
%     [ -1  1  0  0  0  0  0]         [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%     [  0 -1  1  0  0  0  0]         [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%     [  0  0 -1  1  0  0  0]         [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%     [  0  0  0 -1  1  0  0]         [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%     [  0  0  0  0 -1  1  0]         [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%     [  0  0  0  0  0 -1  1]         [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]

L1= eye(m);
L2= zeros(m,m); L2(2:m,1:m-1) = -eye(m-1);
%L3= zeros(m,m); L3(1,1)=-1;
DXB=(L1+L2);

% DYF=[ -1  0  0  0  0 ]   Y*DYF=[ Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14  -Y15 ]
%     [  1 -1  0  0  0 ]         [ Y22-Y21   Y23-Y22   Y24-Y23   Y25-Y24  -Y25 ]
%     [  0  1 -1  0  0 ]         [ Y32-Y31   Y33-Y32   Y34-Y33   Y35-Y34  -Y35 ]
%     [  0  0  1 -1  0 ]         [ Y42-Y41   Y43-Y42   Y44-Y43   Y45-Y44  -Y45 ]
%     [  0  0  0  1 -1 ]         [ Y52-Y51   Y53-Y52   Y54-Y53   Y55-Y54  -Y55 ]
%                                [ Y62-Y61   Y63-Y62   Y64-Y63   Y65-Y64  -Y65 ]
%                                [ Y72-Y71   Y73-Y72   Y74-Y73   Y75-Y74  -Y75 ]

R1= -eye(n);
R2= zeros(n,n); R2(2:n,1:n-1) = eye(n-1);
%R3= zeros(n,n); R3(n,n)=1;
DYF=(R1+R2);

% DYB=[  1 -1  0  0  0 ]   Y*DYB=[ Y11  Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14 ]
%     [  0  1 -1  0  0 ]         [ Y21  Y22-Y21   Y23-Y22   Y24-Y23   Y25-Y24 ]
%     [  0  0  1 -1  0 ]         [ Y32  Y32-Y31   Y33-Y32   Y34-Y33   Y35-Y34 ]
%     [  0  0  0  1 -1 ]         [ Y42  Y42-Y41   Y43-Y42   Y44-Y43   Y45-Y44 ]
%     [  0  0  0  0  1 ]         [ Y52  Y52-Y51   Y53-Y52   Y54-Y53   Y55-Y54 ]
%                                [ Y62  Y62-Y61   Y63-Y62   Y64-Y63   Y65-Y64 ]
%                                [ Y72  Y72-Y71   Y73-Y72   Y74-Y73   Y75-Y74 ]

R1= eye(n);
R2= zeros(n,n); R2(1:n-1,2:n) = -eye(n-1);
%R3= zeros(n,n); R3(1,1)=-1;
DYB=(R1+R2);

