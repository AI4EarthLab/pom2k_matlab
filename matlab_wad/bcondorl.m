function [elf,uaf,vaf,uf,vf,w,ube,ubw] = bcondorl(idx,elf,uaf,vaf,uf,vf,w,ube,ubw,...
                                 im,jm,kb,imm1,jmm1,kbm1,...
                                 fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                                 dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  This is an optional subroutine replacing  bcond and *
% *                using Orlanski's scheme (J. Comp. Phys. 21, 251-269,*
% *                1976), specialized for the seamount problem. To     *
% *                make it work for the seamount problem, I (G.M.)     *
% *                have had to add an extra condition on an "if"       *
% *                statement in the t and s open boundary conditions,  *
% *                which involves the sign of the normal velocity.     *
% *                Thus:                                               *
% *                                                                    *
% *            if(cl==0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k), *
% *                                                                    *
% *                plus 3 others of the same kind.                     *
% *                                                                    *
% **********************************************************************
%
if(idx==1)
    %
    %-----------------------------------------------------------------------
    %
    %     External (2-D) boundary conditions:
    %
    %     In this example the governing boundary conditions are a radiation
    %     condition on uaf(im,j) in the east and an inflow uaf(2,j) in the
    %     west. The tangential velocities are set to zero on both
    %     boundaries. These are only one set of possibilities and may not
    %     represent a choice which yields the most physically realistic
    %     result.
    %
    %     Elevation (in this application, elevation is not a primary
    %     boundary condition):
    %
    for  j=1:jm
        elf(1,j)=elf(2,j);
        elf(im,j)=elf(imm1,j);
    end
    %
    for j=1:jm
        for i=1:im
            elf(i,j)=elf(i,j)*fsm(i,j);
        end
    end
    %
    return
    %
elseif(idx==2)
    %
    %     External (2-D) velocity:
    %
    for j=2:jmm1
        %
        %     West:
        %
        uaf(2,j)=ramp*uabw(j)-sqrt(grav/h(2,j))*(el(2,j)-elw(j));
        uaf(1,j)=uaf(2,j);
        vaf(1,j)=0.e0;
        %
        %     East:
        %
        uaf(im,j)=ramp*uabe(j) ...
            +sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j));
        vaf(im,j)=0.e0;
        %
    end
    %
    for j=1:jm
        for i=1:im
            uaf(i,j)=uaf(i,j)*dum(i,j);
            vaf(i,j)=vaf(i,j)*dvm(i,j);
        end
    end
    %
    return
    %
elseif(idx==3)
    %
    %-----------------------------------------------------------------------
    %
    %     Internal (3-D) boundary conditions:
    %
    %     Eastern and western radiation boundary conditions according to
    %     Orlanski's explicit scheme:
    %
    for k=1:kbm1
        for j=2:jmm1
            %
            %     West:
            %
            denom=(uf(3,j,k)+ub(3,j,k)-2.e0*u(4,j,k));
            if(denom==0.e0)
                denom=0.01e0;
            end
            cl=(ub(3,j,k)-uf(3,j,k))/denom;
            if(cl>1.e0)
                cl=1.e0;
            end
            if(cl<0.e0)
                cl=0.e0;
            end
            uf(2,j,k)=(ub(2,j,k)*(1.e0-cl)+2.e0*cl*u(3,j,k))     ...
                /(1.e0+cl);
            uf(1,j,k)=uf(2,j,k);
            vf(1,j,k)=0.e0;
            %
            %     East:
            %
            denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.e0*u(im-2,j,k));
            if(denom==0.e0)
                denom=0.01e0;
            end
            cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom;
            if(cl>1.e0)
                cl=1.e0;
            end
            if(cl<0.e0)
                cl=0.e0;
            end
            uf(im,j,k)=(ub(im,j,k)*(1.e0-cl)+2.e0*cl*u(im-1,j,k))     ...
                /(1.e0+cl);
            vf(im,j,k)=0.e0;
        end
    end
    %
    for k=1:kbm1
        for j=1:jm
            for i=1:im
                uf(i,j,k)=uf(i,j,k)*dum(i,j);
                vf(i,j,k)=vf(i,j,k)*dvm(i,j);
            end
        end
    end
    %
    return
    %
elseif(idx==4)
    %
    %     Temperature and salinity boundary conditions (using uf and vf:
    %     respectively):
    %
    for k=1:kbm1
        for j=1:jm
            %
            %     West:
            %
            ubw(j,k)=ub(2,j,k);
            denom=(uf(2,j,k)+tb(2,j,k)-2.e0*t(3,j,k));
            if(denom==0.e0)
                denom=0.01e0;
            end
            cl=(tb(2,j,k)-uf(2,j,k))/denom;
            if(cl>1.e0)
                cl=1.e0;
            end
            if(cl<0.e0)
                cl=0.e0;
            end
            uf(1,j,k)=(tb(1,j,k)*(1.e0-cl)+2.e0*cl*t(2,j,k))/(1.e0+cl);
            if(cl==0.e0 && ubw(j,k)>=0.e0)
                uf(1,j,k)=tbw(j,k);
            end
            %
            denom=(vf(2,j,k)+sb(2,j,k)-2.e0*s(3,j,k));
            if(denom==0.e0)
                denom=0.01e0;
            end
            cl=(sb(2,j,k)-vf(2,j,k))/denom;
            if(cl>1.e0)
                cl=1.e0;
            end
            if(cl<0.e0)
                cl=0.e0;
            end
            vf(1,j,k)=(sb(1,j,k)*(1.e0-cl)+2.e0*cl*s(2,j,k))/(1.e0+cl);
            if(cl==0.e0 && ubw(j,k) >= 0.e0)
                vf(1,j,k)=sbw(j,k);
            end
            %
            %     East:
            %
            ube(j,k)=ub(im,j,k);
            denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.e0*t(im-2,j,k));
            if(denom==0.e0) denom=0.01e0;
            end
            cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom;
            if(cl>1.e0) cl=1.e0;
            end
            if(cl<0.e0) cl=0.e0;
            end
            uf(im,j,k)=(tb(im,j,k)*(1.e0-cl)+2.e0*cl*t(im-1,j,k))     ...
                /(1.e0+cl);
            if(cl==0.e0 && ube(j,k)<=0.e0) uf(im,j,k)=tbe(j,k);
            end
            %
            denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.e0*s(im-2,j,k));
            if(denom==0.e0) denom=0.01e0;
            end
            cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom;
            if(cl>1.e0)
                cl=1.e0;
            end
            if(cl<0.e0)
                cl=0.e0;
            end
            vf(im,j,k)=(sb(im,j,k)*(1.e0-cl)+2.e0*cl*s(im-1,j,k))     ...
                /(1.e0+cl);
            if(cl==0.e0 && ube(j,k)<=0.e0)
                vf(im,j,k)=sbe(j,k);
            end
            %
        end
    end
    %
    for k=1:kbm1
        for j=1:jm
            for i=1:im
                uf(i,j,k)=uf(i,j,k)*fsm(i,j);
                vf(i,j,k)=vf(i,j,k)*fsm(i,j);
            end
        end
    end
    
    %
elseif(idx==5)
    %
    %     Vertical velocity boundary conditions:
    %
    for k=1:kbm1
        for j=1:jm
            for i=1:im
                w(i,j,k)=w(i,j,k)*fsm(i,j);
            end
        end
    end
    %
    
    %
elseif(idx==6)
    %
    %     q2 and q2l boundary conditions:
    %
    for k=1:kb
        %
        for j=1:jm
            uf(im,j,k)=1.e-10;
            vf(im,j,k)=1.e-10;
            uf(1,j,k)=1.e-10;
            vf(1,j,k)=1.e-10;
        end
        %
        for j=1:jm
            for i=1:im
                uf(i,j,k)=uf(i,j,k)*fsm(i,j);
                vf(i,j,k)=vf(i,j,k)*fsm(i,j);
            end
        end
    end
    %
    
    %
end
%
end
%
