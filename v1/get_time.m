function [time,ramp] = get_time(iint,time0,lramp,dti,cor)
global pi im jm
    time=dti*iint*1.0/86400.0+time0;
    if(lramp~=0)
        ramp = time/( (2.0*pi)/abs(cor(floor(im/2),floor(jm/2)))/86400.0 );
        if(ramp>1.0)
            ramp=1.0;
        end
    else
        ramp=1.0;
    end
end