function L=OP_AX(m,n)
%compute \overline{f}^x=(f(i,j)+f(i+1,j))/2

% Define additional matrix
% X= [X11 X12 X13 X14 X15]
%    [X21 X22 X23 X24 X25]
%    [X31 X32 X33 X34 X35]
%    [X41 X42 X43 X44 X45]
%    [X51 X52 X53 X54 X55]
%    [X61 X62 X63 X64 X65]
%    [X71 X72 X73 X74 X75]
%
% L=0.5*[1  1  0  0  0  0  0]   L*X=0.5*[ X11+X21  X12+X22  X13+X23  X14+X24  X15+X25 ]
%       [0  1  1  0  0  0  0]           [ X21+X31  X22+X32  X23+X33  X24+X34  X25+X35 ]
%       [0  0  1  1  0  0  0]           [ X31+X41  X32+X42  X33+X43  X34+X44  X35+X45 ]
%       [0  0  0  1  1  0  0]           [ X41+X51  X42+X52  X33+X43  X34+X44  X35+X45 ]
%       [0  0  0  0  1  1  0]           [ X51+X61  X52+X62  X33+X43  X34+X44  X35+X45 ]
%       [0  0  0  0  0  1  1]           [ X61+X71  X62+X72  X33+X43  X34+X44  X35+X45 ]
%       [0  0  0  0  0  1  1]           [ X61+X71  X62+X72  X33+X43  X34+X44  X35+X45 ]

L1= eye(m);
L2= zeros(m,m); L2(1:m-1,2:m) = eye(m-1);
L3= zeros(m,m); L3(m,m-1)=1;
L=0.5*(L1+L2+L3);
end