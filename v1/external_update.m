function [d,d_3d,ua,va,el,uab,vab,elb,egf,vtf,utf,etf,vamax,imax,jmax]=external_update(d,d_3d,uaf,vaf,elf,ua,va,el,uab,vab,elb,iext,egf,vtf,utf,etf)
global smoth isplit kb isp2i ispi h fsm vmaxl;
    
    if(iext==(isplit-2))
                etf=0.25*smoth*elf;
       elseif(iext==(isplit-1))
                etf=etf+0.5*(1.0-0.5*smoth)*elf;
       elseif(iext==isplit)
                etf=(etf+0.5*elf).*fsm;
    end

    [vamax, vapos]=field_max(abs(vaf(:)));
    [imax, jmax]=ind2sub(size(vaf), vapos);    
   
if(vamax<=vmaxl)            
    uab=ua+0.5*smoth*(uab-2.0*ua+uaf);
    vab=va+0.5*smoth*(vab-2.0*va+vaf);
    ua=uaf;
    va=vaf;
    elb=el+0.5*smoth*(elb-2.0*el+elf);
    el=elf;
    d=h+el;
    d_3d=repmat(d,1,1,kb);  %add by hx
    if(iext~=isplit)
        egf=egf+el*ispi;
        utf=utf+2.0* ua .* AXB(d) * isp2i;
        vtf=vtf+2.0* va .* AYB(d) * isp2i;
    end
end
    end