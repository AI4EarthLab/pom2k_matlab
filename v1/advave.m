function  [advua,advva] = advave(aam2d,uab,vab,ua,va,d)
tps = AYB(AXB(d)) .* AXB(AYB(aam2d)) .*  ( DYB(uab)  + DXB(vab)  );
advua = DXB( AXF(AXB(d).*ua) .* AXF(ua)-2.0*d .* aam2d .* DXF(uab) ) ...
        +DYF( ( AXB(AYB(d).*va) .* AYB(ua)-tps ) );
advva = DXF( ( AYB(AXB(d).*ua) .* AXB(va)-tps ) )...
        +DYB( ( AYF(AYB(d).*va) .* AYF(va) - 2.0*d .* aam2d .* DYF(vab) ) );
% if(mode==2)
%     wubot = -AXB(cbc) .* sqrt( uab.^2+AXB(AYF(vab)).^2 ) .* uab;
%     wvbot = -AYB(cbc) .* sqrt( vab.^2+AYB(AXF(uab)).^2 ) .* vab;
%     curv2d = (AYF(va) .* DXC(dy) - AXF(ua) .* DYC(dx))./(dx.*dy);
%     advua = advua-aru  .* AXB( curv2d .* d .* AYF(va) );  
%     advva = advva+arv  .* AYB( curv2d .* d .* AXF(ua) );
% end

return 
