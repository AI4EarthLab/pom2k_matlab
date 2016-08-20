function R=OP_DY(m,n)
% compute \delta_x {f}=(f(i,j)-f(i,j-1))/dx(i,j)

% Define additional matrix
% Y= [Y11 Y12 Y13 Y14 Y15]
%    [Y21 Y22 Y23 Y24 Y25]
%    [Y31 Y32 Y33 Y34 Y35]
%    [Y41 Y42 Y43 Y44 Y45]
%    [Y51 Y52 Y53 Y54 Y55]
%    [Y61 Y62 Y63 Y64 Y65]
%    [Y71 Y72 Y73 Y74 Y75]
%
% R=[  0 -1  0  0  0 ]   R*Y=[ 0  Y12-Y11   Y13-Y12   Y14-Y13   Y15-Y14 ]
%   [  0  1 -1  0  0 ]       [ 0  Y22-Y21   Y23-Y22   Y24-Y23   Y25-Y24 ]
%   [  0  0  1 -1  0 ]       [ 0  Y32-Y31   Y33-Y32   Y34-Y33   Y35-Y34 ]
%   [  0  0  0  1 -1 ]       [ 0  Y42-Y41   Y43-Y42   Y44-Y43   Y45-Y44 ]
%   [  0  0  0  0  1 ]       [ 0  Y52-Y51   Y53-Y52   Y54-Y53   Y55-Y54 ]
%                            [ 0  Y62-Y61   Y63-Y62   Y64-Y63   Y65-Y64 ]
%                            [ 0  Y72-Y71   Y73-Y72   Y74-Y73   Y75-Y74 ]

R1= eye(n);
R2= zeros(n,n); R2(1:n-1,2:n) = -eye(n-1);
R3= zeros(n,n); R3(1,1)=-1;
R=(R1+R2+R3);
end
