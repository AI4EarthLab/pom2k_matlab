
function [elf,uaf,vaf,uf,vf,w] = bcond(idx,elf,uaf,vaf,uf,vf,w,...
                                 im,jm,kb,imm1,jmm1,kbm1,...
                                 fsm,grav,ramp,rfe,h,uabe,ele,el,uabw,rfw,elw,rfn,eln,vabs,rfs,els,...
                                 dum,dvm,hmax,u,v,t,s,tbn,sbn,dti,tbs,sbs,q2,q2l,small,vabn,dx,dy,dt,tbe,sbe,tbw,sbw,zz)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Applies open boundary conditions.                   *
% *                                                                    *
% *                Closed boundary conditions are automatically        *
% *                enabled through specification of the masks, dum,    *
% *                dvm and fsm, in which case the open boundary        *
% *                conditions, included below, will be overwritten.    *
% *                                                                    *
% *                                The C-Grid                          *
% *                                **********                          *
% *                                                                    *
% *                The diagram is for the case where u and v are the   *
% *                primary boundary conditions together with t and     *
% *                s (co-located with el)                              *
% *                                                                    *
% *                All interpolations are centered in space except     *
% *                those at lateral open boundary where an upstream    *
% *                Horizontal locations of e(el), t and s (etc.) are   *
% *                coincident.                                         *
% *                                                                    *
% *                People not acquainted with sigma coordinates have   *
% *                often asked what kind of boundary condition is      *
% *                applied along closed horizontal boundaries.         *
% *                Although the issue is not as important as it might  *
% *                be  for z-level grids, a direct answer is "half-    *
% *                slip" which, of course, is between free slip and    *
% *                non-slip.                                           *
%
%
% East and West end points for the C-grid in POM.
%
%                      west
%
%           v(1,j+1)=0           v(2,j+1) . . . .
%
%     ----<---<----<-----
%     |                 |
%  u(1,j)   el(1,j)   u(2,j)=BC  el(2,j)   u(3,j) . . . .
%             |                   |
%             -----<----<----<-----
%
%           v(1,j)=0              v(2,j) . . . .
%
%                                                    east
%
%                              . . . .  v(im-1,j+1)           v(im,j+1)=0
%
%
%                 . . .  .  u(im-1,j)   el(im-1,j)  u(im,j)=BC  el(im,j)
%                                            |                   |
%                                            ----->----->---->----
%
%                              . . . .   v(im-1,j)             v(im,j)=0
%
%  Notes:
%    1. The suffixes, f  or af, have been deleted.
%    2. All variables NOT designated as boundary condition (=BC) or set to
% zero or obtained from an interior point are calculated points.
%    3. u(1,j) is never used but is obtained from the interior point for
% cosmetic output. Its counterpart, u(im+1,j), fores not exist.
%    4. v=0 at i=1 and i=im are used as open inflow BC s unless specified
% otherwise.
%    5. The south and north extremal points are obtained from the above by
% permuting u to v, v to u, i to j and j to i.


% **********************************************************************

if(idx==1)
    %
    %-----------------------------------------------------------------------
    %
    %     External (2-D) boundary conditions:
    %
    %     In this example, the governing boundary conditions are a radiation
    %     condition on uaf in the east and in the west, and vaf in the north
    %     and south. The tangential velocities are set to zero on both
    %     boundaries. These are only one set of possibilities and may not
    %     represent a choice which yields the most physically realistic
    %     result.
    %
    %     Elevation (in this application, elevation is not a primary
    %     boundary condition):
    %
    
    for j=1:jm
        elf(1,j)=elf(2,j);
        elf(im,j)=elf(imm1,j);
    end
    %
    for i=1:im
        elf(i,1)=elf(i,2);
        elf(i,jm)=elf(i,jmm1);
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
        %     East:
        %
        uaf(im,j)=uabe(j)...
            +rfe*sqrt(grav/h(imm1,j))     ...
            *(el(imm1,j)-ele(j));
        uaf(im,j)=ramp*uaf(im,j);
        vaf(im,j)=0.e0;
        %
        %     West:
        %
        uaf(2,j)=uabw(j)     ...
            -rfw*sqrt(grav/h(2,j))     ...
            *(el(2,j)-elw(j));
        uaf(2,j)=ramp*uaf(2,j);
        uaf(1,j)=uaf(2,j);
        vaf(1,j)=0.0;
        %
    end
    %
    for i=2:imm1
        %
        %     North:
        %
        vaf(i,jm)=vabn(i)     ...
            +rfn*sqrt(grav/h(i,jmm1))     ...
            *(el(i,jmm1)-eln(i));
        vaf(i,jm)=ramp*vaf(i,jm);
        uaf(i,jm)=0.e0;
        %
        %     South:
        %
        vaf(i,2)=vabs(i)     ...
            -rfs*sqrt(grav/h(i,2))     ...
            *(el(i,2)-els(i));
        vaf(i,2)=ramp*vaf(i,2);
        vaf(i,1)=vaf(i,2);
        uaf(i,1)=0.e0;
        %
    end
    %
    %for j=1:jm
    %    for i=1:im
    %        uaf(i,j)=uaf(i,j)*dum(i,j);
    %        vaf(i,j)=vaf(i,j)*dvm(i,j);
    %    end
    %end
    uaf = uaf.*dum;
    vaf = vaf.*dvm;
    %
    return
    %
elseif(idx==3)
    %
    %-----------------------------------------------------------------------
    %
    %     Internal (3-D) boundary conditions:
    %
    %     Velocity (radiation conditions; smoothing is used in the direction
    %     tangential to the boundaries):
    %
    for k=1:kbm1
        for j=2:jmm1
            %
            %     East:
            %
            ga=sqrt(h(im,j)/hmax);
            uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k)     ...
                +.25e0*u(imm1,j+1,k))     ...
                +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k)     ...
                +.25e0*u(im,j+1,k));
            vf(im,j,k)=0.e0;
            %
            %     West:
            %
            ga=sqrt(h(1,j)/hmax);
            uf(2,j,k)=ga*(.25e0*u(3,j-1,k)+.5e0*u(3,j,k)     ...
                +.25e0*u(3,j+1,k))     ...
                +(1.0-ga)*(.25e0*u(2,j-1,k)+.5e0*u(2,j,k)     ...
                +.25e0*u(2,j+1,k));
            uf(1,j,k)=uf(2,j,k);
            vf(1,j,k)=0.e0;
        end
    end
    %
    for k=1:kbm1
        for i=2:imm1
            %
            %     North:
            %
            ga=sqrt(h(i,jm)/hmax);
            vf(i,jm,k)=ga*(0.25*v(i-1,jmm1,k)+0.5*v(i,jmm1,k)     ...
                +0.25*v(i+1,jmm1,k))     ...
                +(1.0-ga)*(0.25*v(i-1,jm,k)+0.5*v(i,jm,k)     ...
                +0.25*v(i+1,jm,k));
            uf(i,jm,k)=0.0;
            %
            %     South:
            %
            ga=sqrt(h(i,1)/hmax);
            vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k)     ...
                +.25e0*v(i+1,3,k))     ...
                +(1.0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k)     ...
                +.25e0*v(i+1,2,k));
            vf(i,1,k)=vf(i,2,k);
            uf(i,1,k)=0.e0;
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
            %     East:
            %
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j));
            if(double(u1<=0.e0))
                uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k));
                vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k));
            else
                uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k));
                vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k));
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(imm1,j));
                    uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1));
                    vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1));
                end
            end
            %
            %     West:
            %
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j));
            if(double(u1>=0.0))
                uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k));
                vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k));
            else
                uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k));
                vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k));
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(2,j));
                    uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1));
                    vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1));
                end
            end
        end
    end
    %
    for k=1:kbm1
        for i=1:im
            %
            %     North:
            %
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1));
            if(double(u1<=0.0))
                uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k));
                vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k));
            else
                uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k));
                vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k));
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(i,jmm1));
                    uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1));
                    vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1));
                end
            end
            %
            %     South:
            %
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2));
            if(double(u1>=0.e0))
                uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k));
         
                vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k));
            else
                uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k));
                vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k));
                if(k~=1 && k~=kbm1)
                    wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti     ...
                        /((zz(k-1)-zz(k+1))*dt(i,2));
                    uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1));
                    vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1));
                end
            end
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
    return
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
    return
    %
elseif(idx==6)
    %
    %     q2 and q2l boundary conditions:
    %
    for k=1:kb
        for j=1:jm
            %
            %     East:
            %
            u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j));
            if(double(u1<=0.e0))
                uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k));
                vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k));
            else
                uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k));
                vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k));
            end
            %
            %     West:
            %
            u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j));
            if(double(u1>=0.e0))
                uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small);
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small);
            else
                uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k));
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k));
            end
        end
    end
    %
    for k=1:kb
        for i=1:im
            %
            %     North:
            %
            u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1));
            if(double(u1<=0.e0))
                uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k));
                vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k));
            else
                uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k));
                vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k));
            end
            %
            %     South:
            %
            u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2));
            if(double(u1>=0.e0))
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small);
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small);
            else
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k));
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k));
            end
        end
    end
    %
    for k=1:kb
        for j=1:jm
            for i=1:im
                uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10;
                vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10;
            end
        end
    end
    %
    return
    %
end
%
