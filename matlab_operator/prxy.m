function  prxy(label,time,a,im,iskp,jm,jskp,scala)

% **********************************************************************
% *                                                                    *
% * FUN%TION    :  Writes a horizontal 2-D field.                      *
% *                                                                    *
% *                label ....... label for output                      *
% *                time ........ time (days)                           *
% *                a(im,jm,kb).. array to be printed                   *
% *                iskp ........ skipping interval for i               *
% *                jskp ........ skipping interval for j               *
% *                scala ....... < 0 for floating point numbers output *
% *                              0 for integer output, divisor for a   *
% *                                based on magnitudes of |a| values   *
% *                              > 0 for integer output, divisor for a *
% *                                given by scala                      *
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
    for j=1:jskp:jm
        for i=1:iskp:im
            amx=max(abs(a(i,j)),amx);
        end
    end
    scale=10.e0^(floor(log10(amx)+100.e0)-103);
end
if(scala>0.e0)
    scale=scala;
end
%
%
fprintf('%s\n',label);
fprintf('Time= %f, multiply all values by %f\n', time,scale);
for ib=1:cols*iskp:im
    %
    ie=ib+(cols-1)*iskp;
    if(ie>im)
        ie=im;
    end
    %
    for j=1:jskp:jm
        jwr=jm+1-j;
        if(scala>=0.e0)
            
            fprintf('%d , %d \n',jwr,floor(a(ib,jwr)/scale+0.5));
            fprintf('%d , %d \n',jwr,floor(a(ie,jwr)/scale+0.5));
            fprintf('%d , %d \n',jwr,floor(a(iskp,jwr)/scale+0.5));
        else
            fprintf('%d , %d \n',jwr,a(ib,jwr));
            fprintf('%d , %d \n',jwr,a(ie,jwr));
            fprintf('%d , %d \n',jwr,a(iskip,jwr));
            
        end
    end
    %
    %
end
return
end
%
