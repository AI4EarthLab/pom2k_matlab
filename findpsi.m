function findpsi(im,jm,imm1,jmm1,time,iskp,jskp,uab,vab,d,dx,dy)
% **********************************************************************
% *                                                                    *
% * ROUTINE NAME:  findpsi                                             *
% *                                                                    *
% * FUN%TION    :  %alculates the stream function, first assuming      *
% *                zero on the southern boundary and then, using the   *
% *                values on the western boundary, the stream function *
% *                is calculated again. If the elevation field is near *
% *                steady state, the two calculations should agree;    *
% *                otherwise not.                                      *
% *                                                                    *
% **********************************************************************
%
%

psi = zeros(im,jm);
%
%     Sweep northward:
%
for j=2:jmm1
    for i=2:im
        psi(i,j+1)=psi(i,j) ...
            +.25e0*uab(i,j)*(d(i,j)+d(i-1,j))     ...
            *(dy(i,j)+dy(i-1,j))
    end
end
%
prxy('Streamfunction, psi from u              ',time,psi,im,iskp,jm,jskp,0.e0);
%
%    Sweep eastward:
%
for j=2:jm
    for i=2:imm1
        psi(i+1,j)=psi(i,j)     ...
            -.25e0*vab(i,j)*(d(i,j)+d(i,j-1))     ...
            *(dx(i,j)+dx(i,j-1))
    end
end
%
prxy('Streamfunction, psi from v              ',time,psi,im,iskp,jm,jskp,0.e0);
end

%
