function [advx,advy,drhox,drhoy,aam] = lateral_viscosity(npg,u,v,dt_3d,aam,aam_init,ub,vb,rho,rmean,ramp,horcon)
global im jm kb dx_3d dy_3d;
    
    [advx,  advy]   = advct (u,v,dt_3d,aam,ub,vb);
    if(npg == 1)  
        [drhox,drhoy]   = baropg(rho, rmean, dt_3d, ramp);
    elseif(npg == 2)
        [drhox,drhoy]   = baropg_mcc(rho, rmean, dt_3d, ramp);
    else
        disp "Error: invalid value for npg";
    end
    
    aam=horcon .* dx_3d .* dy_3d .*sqrt( DXF(u).^2 + DYF(v).^2 +0.5*( DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))) ).^2 );     
    aam(:,1,:)=aam_init;        aam(:,jm,:)=aam_init;
    aam(1,:,:)=aam_init;        aam(im,:,:)=aam_init;
    aam(:,:,kb)=aam_init;

end