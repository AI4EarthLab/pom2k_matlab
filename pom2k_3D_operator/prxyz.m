function prxyz(label,time,a,im,iskp,jm,jskp,kb,ko,nko,scala)

% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Writes horizontal layers of a 3-D field with        *
% *                integers or floating point numbers.                 *
% *                                                                    *
% *                label ....... label for output                      *
% *                time ........ time (days)                           *
% *                a(im,jm,kb).. array to be printed                   *
% *                iskp ........ skipping interval for i               *
% *                jskp ........ skipping interval for j               *
% *                ko .......... 1-D array of k-indices for output     *
% *                nko ......... number of elements in ko              *
% *                scala ....... < 0 for floating point numbers output *
% *                              0 for integer output, divisor for a   *
% *                                based on magnitudes of |a| values   *
% *                              > 0 for integer output, divisor for a *
% *                                given by scala                      *
% *                                                                    *
% *                (NOTE that this combines functions of old prxyz and *
% *                 eprxyz)                                            *
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
if (scala<0.e0) scale = 1.e0;
end
if (scala==0.e0)
    amx=1.e-12;
    for iko=1:nko
        k=ko(iko);
        for j=1:jskp:jm
            for i=1:iskp:im
                amx=max(abs(a(i,j,k)),amx);
            end
        end
    end
    scale=10.e0^(floor(log10(amx)+100.e0)-103);
end
if(scala>0.e0) scale=scala;
end
%
fprintf('%s\n',label);
fprintf('Time= %f, multiply all values by %f\n', time,scale);

for iko=1:nko
    %
    k=ko(iko);
    fprintf('Layer k = %d \n',k);
    %
    for ib=1:cols*iskp:im
        %
        ie=ib+(cols-1)*iskp;
        if(ie>im) ie=im;
            %
            if(scala>=0.e0)
                fprintf('%d        %d            %d\n',ib,ie,iskp);
            else
                fprintf('%d    %d     %d\n',ib,ie,iskp);
            end
            %
            for j=1:jskp:jm
                jwr=jm+1-j;
                if(scala>=0.e0)
                    fprintf('%d , %d \n',jwr,floor(a(ib,jwr,k)/scale+0.5));
                    fprintf('%d , %d \n',jwr,floor(a(ie,jwr,k)/scale+0.5));
                    fprintf('%d , %d \n',jwr,floor(a(iskp,jwr,k)/scale+0.5));
                else
                    fprintf('%d , %d \n',jwr,a(ib,jwr,k));
                    fprintf('%d , %d \n',jwr,a(ie,jwr,k));
                    fprintf('%d , %d \n',jwr,a(iskp,jwr,k));
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
