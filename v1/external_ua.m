function [uaf,advua]=external_ua(h,elb,el,elf,uab,vab,adx2d,drx2d,advua,d,va,e_atmos,wusurf,wubot,iext,aam2d,ua,ramp)
global dte cor grav alpha ispadv uabe uabw rfe rfw ele elw

  if(mod(iext,ispadv)==0)
    advua = DXB( AXF(AXB(d).*ua) .* AXF(ua)-2.0*d .* aam2d .* DXF(uab) ) ...
        +DYF( AXB(AYB(d).*va) .* AYB(ua)- AYB(AXB(d)) .* AXB(AYB(aam2d)) .*  ( DYB(uab)  + DXB(vab)  ) ); 
  end

  uaf=DIVISION((AXB(h+elb) .* uab -2.0* dte .* (adx2d + advua - AXB(cor .* d .* AYF(va)) ...
        + grav.* AXB(d).*( (1.0-2.0*alpha) .* DXB(el) + alpha* (DXB(elb)+ DXB(elf)) + DXB(e_atmos) ) ...
        + drx2d +  (wusurf-wubot) ) ), AXB(h+elf)) ; 
  
  [uaf] = bcond2_ua(uaf,uabe,uabw,rfe,rfw,el,ele,elw,ramp);  
end