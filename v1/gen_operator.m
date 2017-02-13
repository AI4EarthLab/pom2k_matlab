function [OP] = gen_operator(m,n,k)
      OP = Operator();
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



% Suppose:
% X= [X11 X12 X13 X14 X15]      Y= [Y11 Y12 Y13 Y14 Y15]    Z= [Z11 Z12 Z13 Z14 Z15 Z16]
%    [X21 X22 X23 X24 X25]         [Y21 Y22 Y23 Y24 Y25]       [Z21 Z22 Z23 Z24 Z25 Z26]
%    [X31 X32 X33 X34 X35]         [Y31 Y32 Y33 Y34 Y35]       [Z31 Z32 Z33 Z34 Z35 Z36]
%    [X41 X42 X43 X44 X45]         [Y41 Y42 Y43 Y44 Y45]       [Z41 Z42 Z43 Z44 Z45 Z46]
%    [X51 X52 X53 X54 X55]         [Y51 Y52 Y53 Y54 Y55]       [Z51 Z52 Z53 Z54 Z55 Z56]
%    [X61 X62 X63 X64 X65]         [Y61 Y62 Y63 Y64 Y65]       [Z61 Z62 Z63 Z64 Z65 Z66]
%    [X71 X72 X73 X74 X75]         [Y71 Y72 Y73 Y74 Y75]       [Z71 Z72 Z73 Z74 Z75 Z76]


%%%%%%%%%%%%%%%%%%%%%%%% A operator in XY plate %%%%%%%%%%%%%%%%%%
% OP_AXF= 0.5*[1  1  0  0  0  0  0]   OP_AXF*X= 0.5*[ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%             [0  1  1  0  0  0  0]                 [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%             [0  0  1  1  0  0  0]                 [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%             [0  0  0  1  1  0  0]                 [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%             [0  0  0  0  1  1  0]                 [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%             [0  0  0  0  0  1  1]                 [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]
%             [0  0  0  0  0  0  2]                 [    2X71     2X72     2X73     2X74     2X75 ]
L1= eye(m,m); 
L2=[zeros(m,1) eye(m,m-1)];
L3=zeros(m,m); L3(m,m)=1;
OP.OP_AXF=0.5*(L1+L2+L3);

L1= eye(n,n); 
L2=[zeros(n,1) eye(n,n-1)];
L3=zeros(n,n); L3(n,n)=1;
OP.OP_AXF1=0.5*(L1+L2+L3);
%  OP_AXB= 0.5*[2  0  0  0  0  0  0]   OP_AXB*X= 0.5*[    2X11     2X12     2X13     2X14     2X15 ]
%              [1  1  0  0  0  0  0]                 [ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%              [0  1  1  0  0  0  0]                 [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%              [0  0  1  1  0  0  0]                 [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%              [0  0  0  1  1  0  0]                 [ X41+X51  X42+X52  X43+X53  X44+X54  X45+X55 ]
%              [0  0  0  0  1  1  0]                 [ X51+X61  X52+X62  X53+X63  X54+X64  X55+X65 ]
%              [0  0  0  0  0  1  1]                 [ X61+X71  X62+X72  X63+X73  X64+X74  X65+X75 ]
L1= eye(m,m); 
L2=[zeros(1,m); eye(m-1,m)];
L3=zeros(m,m); L3(1,1)=1;
OP.OP_AXB=0.5*(L1+L2+L3);

L1= eye(n,n); 
L2=[zeros(1,n); eye(n-1,n)];
L3=zeros(n,n); L3(1,1)=1;
OP.OP_AXB1=0.5*(L1+L2+L3);
% OP_AYF= 0.5*[1  0  0  0  0]        Y*OP_AYF= 0.5*[ Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15  2Y15 ]
%             [1  1  0  0  0]                      [ Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25  2Y25 ]
%             [0  1  1  0  0]                      [ Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35  2Y35 ]
%             [0  0  1  1  0]                      [ Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45  2Y45 ]
%             [0  0  0  1  2]                      [ Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55  2Y55 ]
%                                                  [ Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65  2Y65 ]
%                                                  [ Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75  2Y75 ]
R1= eye(n,n); 
R2=[zeros(1,n); eye(n-1,n)];
R3=zeros(n,n); R3(n,n)=1;
OP.OP_AYF=0.5*(R1+R2+R3);



% OP_AYB= 0.5*[2  1  0  0  0]        Y*OP_AYB= 0.5*[ 2Y11  Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15 ]
%             [0  1  1  0  0]                      [ 2Y21  Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25 ]
%             [0  0  1  1  0]                      [ 2Y31  Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35 ]
%             [0  0  0  1  1]                      [ 2Y41  Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45 ]
%             [0  0  0  0  1]                      [ 2Y51  Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55 ]
%                                                  [ 2Y61  Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65 ]
%                                                  [ 2Y71  Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75 ]
R1= eye(n,n); 
R2=[zeros(n,1), eye(n,n-1)];
R3=zeros(n,n); R3(1,1)=1;
OP.OP_AYB=0.5*(R1+R2+R3);

%%%%%%%%%%%%%%%%%%%%%%%% A operator in XZ plate %%%%%%%%%%%%%%%%%%
R1= eye(k,k); 
R2=[zeros(1,k); eye(k-1,k)];
R3=zeros(k,k); R3(k,k)=1;
OP.OP_AZF=0.5*(R1+R2+R3);

R1= eye(k,k); 
R2=[zeros(k,1), eye(k,k-1)];
R3=zeros(k,k); R3(1,1)=1;
OP.OP_AZB=0.5*(R1+R2+R3);


R1 = [zeros(k,1), eye(k, k-1)];
R2 = [zeros(1,k); -eye(k-1, k)];
OP.OP_Rk = R1 + R2;

%%%%%%%%%%%%%%%%%%%%%%%% D operator in XY plate %%%%%%%%%%%%%%%%%%
% OP_DXF1_XY=[ -1  1  0  0  0  0  0]    OP_DXF1_XY*X =[ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%            [  0 -1  1  0  0  0  0]                  [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%            [  0  0 -1  1  0  0  0]                  [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%            [  0  0  0 -1  1  0  0]                  [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%            [  0  0  0  0 -1  1  0]                  [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%            [  0  0  0  0  0 -1  1]                  [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
%            [  0  0  0  0  0  0  0]                  [       0        0         0         0         0 ]
L1= -eye(m,m); 
L2=[zeros(m,1) eye(m,m-1)];
L3=zeros(m,m); L3(m,m)=1;
OP.OP_DXF=(L1+L2+L3);

% OP_DXB1_XY=[  0  0  0  0  0  0  0]    OP_DXB1_XY*X =[       0        0         0         0         0 ]
%            [ -1  1  0  0  0  0  0]                  [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%            [  0 -1  1  0  0  0  0]                  [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%            [  0  0 -1  1  0  0  0]                  [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%            [  0  0  0 -1  1  0  0]                  [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%            [  0  0  0  0 -1  1  0]                  [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%            [  0  0  0  0  0 -1  1]                  [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]
L1=eye(m,m); 
L2=[zeros(1,m); -eye(m-1,m)];
L3=zeros(m,m); L3(1,1)=-1;
OP.OP_DXB=(L1+L2+L3);

% OP_DYF1_XY=[ -1  0  0  0  0 ]         Y*OP_DYF1_XY =[ Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14  0 ]
%            [  1 -1  0  0  0 ]                       [ Y22-Y21   Y23-Y22   Y24-Y23   Y25-Y24  0 ]
%            [  0  1 -1  0  0 ]                        [ Y22-Y21   Y33-Y32   Y34-Y33   Y35-Y34  0 ]
%            [  0  0  1 -1  0 ]                       [ Y22-Y21   Y43-Y42   Y44-Y43   Y45-Y44  0 ]
%            [  0  0  0  1  0 ]                       [ Y22-Y21   Y53-Y52   Y54-Y53   Y55-Y54  0 ]
%                                                     [ Y22-Y21   Y63-Y62   Y64-Y63   Y65-Y64  0 ]
%                                                     [ Y22-Y21   Y73-Y72   Y74-Y73   Y75-Y74  0 ]
R1= -eye(n,n); 
R2=[zeros(1,n); eye(n-1,n)];
R3=zeros(n,n); R3(n,n)=1;
OP.OP_DYF=(R1+R2+R3);

% OP_DYB1_XY=[  0 -1  0  0  0 ]         Y*OP_DYB1_XY =[ 0  Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14 ]
%            [  0  1 -1  0  0 ]                       [ 0  Y22-Y21   Y23-Y22   Y24-Y23   Y15-Y14 ]
%            [  0  0  1 -1  0 ]                       [ 0  Y32-Y31   Y33-Y32   Y34-Y33   Y35-Y34 ]
%            [  0  0  0  1 -1 ]                       [ 0  Y42-Y41   Y43-Y42   Y44-Y43   Y45-Y44 ]
%            [  0  0  0  0  1 ]                       [ 0  Y52-Y51   Y53-Y52   Y54-Y53   Y55-Y54 ]
%                                                     [ 0  Y62-Y61   Y63-Y62   Y64-Y63   Y65-Y64 ]
%                                                     [ 0  Y72-Y71   Y73-Y72   Y74-Y73   Y75-Y74 ]
R1= eye(n,n); 
R2=[zeros(n,1), -eye(n,n-1)];
R3=zeros(n,n); R3(1,1)=-1;
OP.OP_DYB=(R1+R2+R3);


R1= -eye(k,k); 
R2=[zeros(1,k); eye(k-1,k)];
R3=zeros(k,k); R3(k,k)=1;
OP.OP_DZF=(R1+R2+R3);

R1= eye(k,k); 
R2=[zeros(k,1), -eye(k,k-1)];
R3=zeros(k,k); R3(1,1)=-1;
OP.OP_DZB=(R1+R2+R3);
%%%%%%%%%%%%%%%%%%%%%%%% OP_L and OP_R operator %%%%%%%%%%%%%%%%%%
L=eye(m,m); L(1,1)=0;L(m,m)=0;
R=eye(n,n); R(1,1)=0;R(n,n)=0;
OP.OP_L_XY = L;
OP.OP_R_XY = R;

L=eye(m,m); L(1,1)=0;L(m,m)=0;
R=eye(k,k); R(1,1)=0;R(k,k)=0;
OP.OP_L_XZ = L;
OP.OP_R_XZ = R;

L=eye(n,n); L(1,1)=0;L(n,n)=0;
R=eye(k,k); R(1,1)=0;R(k,k)=0;
OP.OP_L_YZ = L;
OP.OP_R_YZ = R;

%%%%%%%%%%%%%%%%%%%%%%%% OP_SUM operator %%%%%%%%%%%%%%%%%%
% Define sum operator along with X and Y directions. 
%   R1=[1 1 1 1 1 0]     R2=[0 1 1 1 1 1]
%      [0 1 1 1 1 0]        [0 0 1 1 1 1] 
%      [0 0 1 1 1 0]        [0 0 0 1 1 1]
%      [0 0 0 1 1 0]        [0 0 0 0 1 1]
%      [0 0 0 0 1 0]        [0 0 0 0 0 1]
%      [0 0 0 0 0 0]        [0 0 0 0 0 0]

R1=triu(ones(k,k)); R1(:,k)=0;
OP.OP_SUMZ1= R1;

R2=triu(ones(k,k))-eye(k);
OP.OP_SUMZ2= R2;

end

