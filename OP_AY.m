function R=OP_AY(m,n)
%compute \overline{f}^y=(f(i,j)+f(i,j+1))/2
 
% Define additional matrix
% Y= [Y11 Y12 Y13 Y14 Y15]
%    [Y21 Y22 Y23 Y24 Y25]
%    [Y31 Y32 Y33 Y34 Y35]
%    [Y41 Y42 Y43 Y44 Y45]
%    [Y51 Y52 Y53 Y54 Y55]
%    [Y61 Y62 Y63 Y64 Y65]
%    [Y71 Y72 Y73 Y74 Y75]
%
% R=0.5*[1  0  0  0  0]   Y*R=0.5*[ Y11+Y12  Y12+Y13  Y13+Y14  Y14+Y15  Y14+Y15 ]
%       [1  1  0  0  0]           [ Y21+Y22  Y22+Y23  Y23+Y24  Y24+Y25  Y24+Y25 ]
%       [0  1  1  0  0]           [ Y31+Y32  Y32+Y33  Y33+Y34  Y34+Y35  Y34+Y35 ]
%       [0  0  1  1  1]           [ Y41+Y42  Y42+Y43  Y43+Y44  Y44+Y45  Y44+Y45 ]
%       [0  0  0  1  1]           [ Y51+Y52  Y52+Y53  Y53+Y54  Y54+Y55  Y54+Y55 ]
%                                 [ Y61+Y62  Y62+Y63  Y63+Y64  Y64+Y65  Y64+Y65 ]
%                                 [ Y71+Y72  Y72+Y73  Y73+Y74  Y74+Y75  Y74+Y75 ]
R1= eye(n);
R2= zeros(n,n); R2(2:n,1:n-1)= eye(n-1);
R3= zeros(n,n); R3(n-1,n)=1;
R = 0.5*(R1+R2+R3);
end
