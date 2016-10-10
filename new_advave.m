function  [curv2d,wubot,wvbot,advua,advva] = new_advave(curv2d,wubot,wvbot,mode,aam2d,uab,vab,ua,va,cbc,d)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates horizontal advection and diffusion.      *
% *                                                                    *
% **********************************************************************
load('grid.mat');load('xyz.mat');


tps = zeros(im,jm);
tps = AYB1_XY(AXB1_XY(d)) .* AXB1_XY(AYB1_XY(aam2d)) .* ( DYB1_XY(uab)...
     .*DIVISION( 1.0,AYB1_XY(AXB1_XY(dy)) ) + DXB1_XY(vab) .* DIVISION( 1.0,AYB1_XY(AXB1_XY(dx)) ) );

advua = zeros(im,jm);
advva=zeros(im,jm);
advua = DXB2_XY(( AXF1_XY(AXB1_XY(d).*ua) .* AXF2_XY(ua)-2.0*d .* aam2d .* DXF2_XY(uab)./dx ) .* dy ) ...
        +DYF2_XY( ( AXB1_XY(AYB1_XY(d).*va) .* AYB1_XY(ua)-tps ) .* AYB1_XY(AXB1_XY(dx)) );
advva = DXF2_XY( ( AYB1_XY(AXB1_XY(d).*ua) .* AXB1_XY(va)-tps ) .* AYB1_XY(AXB1_XY(dy)) )...
        +DYB2_XY( ( AYF1_XY(AYB1_XY(d).*va) .* AYF1_XY(va) - 2.0*d .* aam2d .* DYF2_XY(vab)./dy ) .* dx );

%%%%%%%%%%%%%%%%%%% 
%wubot = zeros(im,jm);
%wvbot = zeros(im,jm);

if(mode==2)
    wubot = -AXB2_XY(cbc) .* sqrt( uab.^2+AXB2_XY(AYF2_XY(vab)).^2 ) .* uab;
    wvbot = -AYB2_XY(cbc) .* sqrt( vab.^2+AYB2_XY(AXF2_XY(uab)).^2 ) .* vab;
    curv2d = (AYF2_XY(va) .* DXC_XY(dy) - AXF2_XY(ua) .* DYC_XY(dx))./(dx.*dy);
    advua = advua-aru .* AXB2_XY( curv2d .* d .* AYF2_XY(va) );  
    advva = advva+arv .* AYB2_XY( curv2d .* d .* AXF2_XY(ua) );
end

return 