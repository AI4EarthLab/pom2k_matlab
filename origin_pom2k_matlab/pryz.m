function pryz(label,time,a,im,jm,jskp,kb,io,nio,scala,dt,zz)

% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
% *                y- or j-direction.                                  *
% *                                                                    *
% *                label ....... label for output                      *
% *                time ........ time (days)                           *
% *                a(im,jm,kb).. array to be printed                   *
% *                jskp ........ skipping interval for j               *
% *                io .......... 1-D array of i-indices for output     *
% *                nio ......... number of elements in io              *
% *                scala ....... < 0 for floating point numbers output *
% *                              0 for integer output, divisor for a   *
% *                                based on magnitudes of |a| values   *
% *                              > 0 for integer output, divisor for a *
% *                                given by scala                      *
% *                dt(im,jm) ... total depth                           *
% *                zz(kb) ...... sigma coordinate                      *
% *                                                                    *
% **********************************************************************
%
%
if(scala>=0.e0)
    cols=24;
else
    cols=12;
end
%
if (scala<0.e0)
    scale = 1.e0;
end
if (scala==0.e0)
    amx=1.e-12
    for  k=1:kb
        for j=1:jskp:jm
            for iio=1:nio
                i=io(iio);
                amx=max(abs(a(i,j,k)),amx);
            end
        end
    end
    scale=10.e0^(floor(log10(amx)+100.e0)-103);
end
if(scala>0.e0)
    scale=scala;
end
%

fprintf('%s\n',label);
fprintf('Time= %f, multiply all values by %f\n', time,scale);
%
for iio=1:nio
    %
    i=io(iio);
    %
    fprintf('Section i = %d \n',i);
    for jb=1:cols*jskp:jm
        %
        je=jb+(cols-1)*jskp;
        if(je>jm)
            je=jm;
        end
        %
        if(scala>=0.e0)
            fprintf('%20d%20d%20d',jb,je,jskp);
        else
            fprintf('%10d%10d%10d',jb,je,jskp);
        end
        
        fprintf('%10d \n',floor(dt(i,jb)/scale+0.5));
        fprintf('%10d \n',floor(dt(i,je)/scale+0.5));
        fprintf('%10d \n',floor(dt(i,jskp)/scale+0.5));
        
        for k=1:kb
            if(scala>=0.e0)
                fprintf('%d %d %d \n',k,zz(k),floor(a(i,jb,k)/scale+0.5));
                fprintf('%d %d %d \n',k,zz(k),floor(a(i,je,k)/scale+0.5));
                fprintf('%d %d %d \n',k,zz(k),floor(a(i,jskp,k)/scale+0.5));
                
            else
                fprintf('%d %d %d \n',k,zz(k),a(i,jb,k));
                fprintf('%d %d %d \n',k,zz(k),a(i,je,k));
                fprintf('%d %d %d \n',k,zz(k),a(i,jb,jskp));
            end
        end
        %
        fprintf('//');
        %
    end
    %
end
%
return
%
end
%
