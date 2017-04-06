function [drhox,drhoy]=baropg_mcc(rho,rmean,d_3d,dt_3d,ramp)
% calculate  baroclinic pressure gradient 4th order correction terms, following McCalpin
global im kb dum dvm dx_3d dy_3d dz_3d grav zz_3d gs;
dum_3d=repmat(dum,1,1,kb);    dvm_3d=repmat(dvm,1,1,kb);

rho=rho-rmean;

    tmpdrho=  DXB( rho) .* dum_3d;      tmprhou=  AXB(rho) .*dum_3d;
    tmpddx=   DXB( d_3d).* dum_3d;      tmpd4=    AXB(d_3d).*dum_3d;
    tmpdrho=create_field(double(tmpdrho),gs,3);     

    drho=tmpdrho -( DXB ( DXF ( DXB( double(rho) ) .* double(dum_3d) ) )./ AXB(dx_3d)-2.0.*DXB( rho ) .* (1-dum_3d) ) ./24.e0;
    rhou=tmprhou-DXB(2.0 .* AXF(DXB( double(rho) ) .* double(dum_3d) ) )./16.e0;
    ddx =tmpddx -( DXB ( DXF  ( DXB( double(d_3d)) .* double(dum_3d) ) )./ AXB(dx_3d)-2.0.*DXB( d_3d) .* (1-dum_3d) ) ./24.e0;
    d4  =tmpd4-  DXB(2.0 .* AXF(DXB( double(d_3d)) .* double(dum_3d) )) ./16.e0;
    drho(1:2,:,:)   =tmpdrho(1:2,:,:);       rhou(1:2,:,:)  =tmprhou(1:2,:,:);
    ddx(1:2,:,:)    =tmpddx(1:2,:,:);        d4(1:2,:,:)    =tmpd4(1:2,:,:);
        
drhox=AXB(dt_3d) .* ramp .* grav .* SUMZ(-DZB(zz_3d).* d4 .* AZB(drho) +  AZB(zz_3d) .* ddx .* DZB(rhou).*AZB(dz_3d)).* dum_3d;
drhox(im,:,:)=0.e0;
    tmpdrho=   DYB( rho) .* dvm_3d;      tmprhou=   AYB(rho).* dvm_3d;
    tmpddx=    DYB( d_3d).*dvm_3d;       tmpd4=     AYB(d_3d).*dvm_3d;      
    tmpdrho=create_field(double(tmpdrho),gs,3);
        
    drho=tmpdrho-(DYB ( DYF ( DYB( double(rho) ) .* double(dvm_3d) ) )./ AYB(dy_3d) - 2.0.*DYB( rho)  .*(1-dvm_3d) )./24.e0;
    rhou=tmprhou-DYB(2.0 .* AYF(DYB( double(rho) )  .* double(dvm_3d) ) )./16.e0;
    ddx =tmpddx -(DYB ( DYF ( DYB( double(d_3d)) .* double(dvm_3d) ) )./ AYB(dy_3d) - 2.0.*DYB( d_3d) .*(1-dvm_3d) )./24.e0;
    d4  =tmpd4  -DYB(2.0 .* AYF(DYB( double(d_3d) ) .* double(dvm_3d) ) )./16.e0;
    drho(:,1:2,:)   =tmpdrho(:,1:2,:);       rhou(:,1:2,:)  =tmprhou(:,1:2,:);
    ddx(:,1:2,:)    =tmpddx(:,1:2,:);        d4(:,1:2,:)    =tmpd4(:,1:2,:);        
                
drhoy=AYB(dt_3d) .* ramp .* grav .* SUMZ(-DZB(zz_3d).* d4 .* AZB(drho) +  AZB(zz_3d) .* ddx .* DZB(rhou).*AZB(dz_3d)).* dvm_3d;
drhoy(im,:,:)=0.e0;              
end