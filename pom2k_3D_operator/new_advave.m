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

tps = AYB1(AXB1(d)) .* AXB1(AYB1(aam2d)) ...
       .*  ( DYB1(uab) .*DIVISION( 1.0,AYB1(AXB1(dy)) ) + DXB1(vab) .* DIVISION( 1.0,AYB1(AXB1(dx)) ) );

advua = DXB2(( AXF1(AXB1(d).*ua) .* AXF2(ua)-2.0*d .* aam2d .* DXF2(uab)./dx ) .* dy ) ...
        +DYF2( ( AXB1(AYB1(d).*va) .* AYB1(ua)-tps ) .* AYB1(AXB1(dx)) );
advva = DXF2( ( AYB1(AXB1(d).*ua) .* AXB1(va)-tps ) .* AYB1(AXB1(dy)) )...
        +DYB2( ( AYF1(AYB1(d).*va) .* AYF1(va) - 2.0*d .* aam2d .* DYF2(vab)./dy ) .* dx );

%%%%%%%%%%%%%%%%%%% 
%wubot = zeros(im,jm);
%wvbot = zeros(im,jm);

if(mode==2)
    wubot = -AXB2(cbc) .* sqrt( uab.^2+AXB2(AYF2(vab)).^2 ) .* uab;
    wvbot = -AYB2(cbc) .* sqrt( vab.^2+AYB2(AXF2(uab)).^2 ) .* vab;
    curv2d = (AYF2(va) .* DXC(dy) - AXF2(ua) .* DYC(dx))./(dx.*dy);
    advua = advua-aru  .* AXB2( curv2d .* d .* AYF2(va) );  
    advva = advva+arv  .* AYB2( curv2d .* d .* AXF2(ua) );
end

return 
