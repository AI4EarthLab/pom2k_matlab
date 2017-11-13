function prxz(label,time,a,im,iskp,jm,kb,jo,njo,scala,dt,zz)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Writes vertical section of a 3-D field, in the      *
% *                x- or i-direction .                                 *
% *                                                                    *
% *                label ....... label for output                      *
% *                time ........ time (days)                           *
% *                a(im,jm,kb).. array to be printed                   *
% *                iskp ........ skipping interval for i               *
% *                jo .......... 1-D array of j-indices for output     *
% *                njo ......... number of elements in jo              *
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
    amx=1.e-12;
    for  k=1:kb
        for ijo=1:njo
            j=jo(ijo);
            for i=1:iskp:im
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

for ijo=1:njo
    j=jo(ijo);
    fprintf('Section j = %d \n',j);
    for ib=1:cols*iskp:im
        ie=ib+(cols-1)*iskp;
        if(ie>im) ie=im;
        end
        
        fprintf('i = %d , d= %d , z or zz',ib,floor(dt(ib,j)+0.5));
        fprintf('i = %d , d= %d , z or zz',ie,floor(dt(ie,j)+0.5));
        fprintf('i = %d , d= %d , z or zz',iskp,floor(dt(iskp,j)+0.5));
        
        for k=1:kb
            if(scala>=0.e0)
                fprintf('%d , %f  %d\n',k,zz(k),floor(a(ib,j,k)/scale+0.5));
                fprintf('%d , %f  %d\n',k,zz(k),floor(a(ie,j,k)/scale+0.5));
                fprintf('%d , %f  %d\n',k,zz(k),floor(a(iskp,j,k)/scale+0.5));
            else
                fprintf('%d , %f  %f\n',k,zz(k),a(ib,j,k));
                fprintf('%d , %f  %f\n',k,zz(k),a(ie,j,k));
                fprintf('%d , %f  %f\n',k,zz(k),a(ik,j,k));
            end
        end
        %
    end
    %
end
%
return
%
end
%
