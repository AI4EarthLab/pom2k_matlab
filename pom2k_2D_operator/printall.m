function printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s,rho,aam,km,kb,mode,dt,zz,z)
% **********************************************************************
% *                                                                    *
% *                         POM2K SOURCE CODE                          *
% *                                                                    *
% * ROUTINE NAME:  printall                                            *
% *                                                                    *
% * FUNCTION    :  Prints a set of outputs to device 6                 *
% *                                                                    *
% *                Edit as approriate.                                 *
% *                                                                    *
% **********************************************************************
%
%
%     2-D horizontal fields:
%
prxy('Depth-averaged u, uab                   ',time,uab,im,iskp,jm,jskp,0.e0);
%
prxy('Depth-averaged v, vab                   ',time,vab,im,iskp,jm,jskp,0.e0);
%
prxy('Surface elevation, elb                  ',time,elb,im,iskp,jm,jskp,0.e0);
%
%
%     Calculate and print streamfunction:
%
findpsi(im,jm,imm1,jmm1,time,iskp,jskp,uab,vab,d,dx,dy);
%
if(mode~=2)
    %
    %     2-D horizontal sections of 3-D fields:
    %
    %     Set levels for output:
    %
    ko(1)=1;
    ko(2)=floor(kb/2);
    ko(3)=kb-1;
    %
    prxyz('x-velocity, u                           ',time,u,im,iskp,jm,jskp,kb,ko,3,0.e0 );
    %
    prxyz('y-velocity, v                           ',time,v,im,iskp,jm,jskp,kb,ko,3,0.e0 );
    %
    ko(1)=2;
    prxyz('z-velocity, w                           ',time,w,im,iskp,jm,jskp,kb,ko,3,0.e0 );
    ko(1)=1;
    %
    prxyz('Potential temperature, t                ',time,t,im,iskp,jm,jskp,kb,ko,3,1.e-2);
    %
    prxyz('Salinity, s                              ',time,s,im,iskp,jm,jskp,kb,ko,3,1.e-2);
    %
    prxyz('(density-1000)/rhoref, rho              ',time,rho,im,iskp,jm,jskp,kb,ko,3,1.e-5);
    %
    %
    prxyz('Horizontal kinematic viscosity, aam     ',time,aam,im,iskp,jm,jskp,kb,ko,3,0.e0 );
    %
    prxyz('Vertical kinematic viscosity, km        ',time,km,im,iskp,jm,jskp,kb,ko,3,0.e0 );
    %
    %
    %     Set sections for output:
    %
    jo(1)=1;
    jo(2)=floor(jm/2);
    jo(3)=jm-1;
    
    prxz('x-velocity, u                           ',time,u,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz);
    %
    prxz('y-velocity, v                           ',time,v,im,iskp,jm,kb,jo,3,0.e0 ,dt,zz);
    %
    prxz('z-velocity, w                           ',time,w,im,iskp,jm,kb,jo,3,0.e0 ,dt,z );
    %
    prxz('Potential temperature, t                ',time,t,im,iskp,jm,kb,jo,3,1.e-2,dt,zz);
    %
    prxz('Salinity, s                             ',time,s,im,iskp,jm,kb,jo,3,1.e-2,dt,zz);
    %
    prxz('(density-1000)/rhoref, rho              ',time,rho,im,iskp,jm,kb,jo,3,1.e-5,dt,zz);
    %
    %
    %     Set se%tions for output:
    %
    io(1)=1;
    io(2)=floor(im/2);
    io(3)=im-1;
    %
    pryz('x-velocity, u                           ', time,u,im,jm,jskp,kb,io,3,0.e0 ,dt,zz);
    %
    pryz('y-velocity, v                           ', time,v,im,jm,jskp,kb,io,3,0.e0 ,dt,zz);
    %
    pryz('z-velocity, w                           ', time,w,im,jm,jskp,kb,io,3,0.e0 ,dt,zz);
    %
    pryz('Potential temperature, t                ', time,t,im,jm,jskp,kb,io,3,1.e-2,dt,zz);
    %
end
%
