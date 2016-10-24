function  [curv2d,wubot,wvbot,advua,advva] = new_advave(curv2d,wubot,wvbot,mode,aam2d,uab,vab,ua,va,cbc,d)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates horizontal advection and diffusion.      *
% *                                                                    *
% **********************************************************************
load('grid.mat');

advua = zeros(im,jm);
advva = zeros(im,jm);
tps = zeros(im,jm);

tps = AYB(AXB(d)) .* AXB(AYB(aam2d)) ...
       .*  ( DYB(uab) .*DIVISION( 1.0,AYB(AXB(dy)) ) + DXB(vab) .* DIVISION( 1.0,AYB(AXB(dx)) ) );

advua = DXB(( AXF(AXB(d).*ua) .* AXF(ua)-2.0*d .* aam2d .* DXF(uab)./dx ) .* dy ) ...
        +DYF( ( AXB(AYB(d).*va) .* AYB(ua)-tps ) .* AYB(AXB(dx)) );
advva = DXF( ( AYB(AXB(d).*ua) .* AXB(va)-tps ) .* AYB(AXB(dy)) )...
        +DYB( ( AYF(AYB(d).*va) .* AYF(va) - 2.0*d .* aam2d .* DYF(vab)./dy ) .* dx );

%%%%%%%%%%%%%%%%%%% 
%wubot = zeros(im,jm);
%wvbot = zeros(im,jm);

if(mode==2)
    wubot = -AXB(cbc) .* sqrt( uab.^2+AXB(AYF(vab)).^2 ) .* uab;
    wvbot = -AYB(cbc) .* sqrt( vab.^2+AYB(AXF(uab)).^2 ) .* vab;
    curv2d = (AYF(va) .* DXC(dy) - AXF(ua) .* DYC(dx))./(dx.*dy);
    advua = advua-aru  .* AXB( curv2d .* d .* AYF(va) );  
    advva = advva+arv  .* AYB( curv2d .* d .* AXF(ua) );
end

return 
