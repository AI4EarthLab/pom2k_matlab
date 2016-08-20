function L=OP_DX(m,n)
% compute \delta_x {f}=(f(i,j)-f(i-1,j))/2

% Define additional matrix
% X= [X11 X12 X13 X14 X15]
%    [X21 X22 X23 X24 X25]
%    [X31 X32 X33 X34 X35]
%    [X41 X42 X43 X44 X45]
%    [X51 X52 X53 X54 X55]
%    [X61 X62 X63 X64 X65]
%    [X71 X72 X73 X74 X75]
%
% L=[  0  0  0  0  0  0  0]   L*X=[       0         0         0         0         0 ]
%   [ -1  1  0  0  0  0  0]       [ X21-X11  X22-X12   X23-X13   X24-X14   X25-X15 ]
%   [  0 -1  1  0  0  0  0]       [ X31-X21  X32-X22   X33-X23   X34-X24   X35-X25 ]
%   [  0  0 -1  1  0  0  0]       [ X41-X31  X42-X32   X43-X33   X44-X34   X45-X35 ]
%   [  0  0  0 -1  1  0  0]       [ X51-X41  X52-X42   X53-X43   X54-X44   X55-X45 ]
%   [  0  0  0  0 -1  1  0]       [ X61-X51  X62-X52   X63-X53   X64-X54   X65-X55 ]
%   [  0  0  0  0  0 -1  1]       [ X71-X61  X72-X62   X73-X63   X74-X64   X75-X65 ]

L1= eye(m);
L2= zeros(m,m); L2(2:m,1:m-1) = -eye(m-1);
L3= zeros(m,m); L3(1,1)=-1;
L=(L1+L2+L3);
end
