function [vaf,advva] = external_va(h,elb,el,elf,uab,vab,ady2d,dry2d,advva,d,ua,e_atmos,wvsurf,wvbot,iext,aam2d,va,ramp)
global dte cor grav alpha ispadv vabn vabs rfn rfs eln els 

    if(mod(iext,ispadv)==0)
      advva = DXF( AYB(AXB(d).*ua) .* AXB(va)- AYB(AXB(d)) .* AXB(AYB(aam2d)) .* ( DYB(uab) + DXB(vab) ) )...
        +DYB( ( AYF(AYB(d).*va) .* AYF(va) - 2.0*d .* aam2d .* DYF(vab) ) );
    end

vaf=   DIVISION((AYB(h+elb) .* vab -2.0* dte .* (ady2d + advva + AYB(cor .* d .* AXF(ua)) ...
    + grav .* AYB(d).*( (1.0-2.0*alpha) .* DYB(el) + alpha* (DYB(elb)+ DYB(elf)) + DYB(e_atmos) ) ...
    + dry2d + (wvsurf-wvbot) )), AYB(h+elf)) ;

[vaf] = bcond2_va(vaf,vabn,vabs,rfn,rfs,el,eln,els,ramp);
end