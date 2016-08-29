
function [curv2d,advua,advva,fluxua,fluxva,wubot,wvbot,tps] = advave(curv2d,advua,advva,fluxua,fluxva,wubot,wvbot,tps,...
                                                                     mode,im,jm,imm1,jmm1,aam2d,...
                                                                     uab,vab,dx,dy,ua,va,cbc,aru,arv,d)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Calculates horizontal advection and diffusion.      *
% *                                                                    *
% **********************************************************************
%
%
%
%
%     u-advection and diffusion:
%
%     Advective fluxes:
%



advua=zeros(im,jm);



%
for j=2:jm
    for i=2:imm1
        fluxua(i,j)=0.125*((d(i+1,j)+d(i,j))*ua(i+1,j) ...
            +(d(i,j)+d(i-1,j))*ua(i,j)) ...
            *(ua(i+1,j)+ua(i,j));
    end
end
%
for j=2:jm
    for i=2:im
        fluxva(i,j)=0.125*((d(i,j)+d(i,j-1))*va(i,j)    ...
            +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))...
            *(ua(i,j)+ua(i,j-1));
    end
end
%
%     Add viscous fluxes:
%
for j=2:jm
    for i=2:imm1
        fluxua(i,j)=fluxua(i,j)     ...
            -d(i,j)*2*aam2d(i,j)*(uab(i+1,j)-uab(i,j))     ...
            /dx(i,j);
    end
end
%
for j=2:jm
    for i=2:im
        tps(i,j)=0.25*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))...
            *(aam2d(i,j)+aam2d(i,j-1)...
            +aam2d(i-1,j)+aam2d(i-1,j-1))...
            *((uab(i,j)-uab(i,j-1))...
            /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))...
            +(vab(i,j)-vab(i-1,j))...
            /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)));
        fluxua(i,j)=fluxua(i,j)*dy(i,j);
        fluxva(i,j)=(fluxva(i,j)-tps(i,j))*0.25     ...
            *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1));
    end
end
%
for j=2:jmm1
    for i=2:imm1
        advua(i,j)=fluxua(i,j)-fluxua(i-1,j)     ...
            +fluxva(i,j+1)-fluxva(i,j);
    end
end
%
%     v-advection and diffusion:

advva=zeros(im,jm);
%
%     Advective fluxes:
%
for j=2:jm
    for i=2:im
        fluxua(i,j)=0.125*((d(i,j)+d(i-1,j))*ua(i,j)     ...
            +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))     ...
            *(va(i-1,j)+va(i,j));
    end
end
%
for j=2:jmm1
    for i=2:im
        fluxva(i,j)=0.125*((d(i,j+1)+d(i,j))*va(i,j+1)     ...
            +(d(i,j)+d(i,j-1))*va(i,j))     ...
            *(va(i,j+1)+va(i,j));
    end
end
%
%     Add viscous fluxes:
%
for j=2:jmm1
    for i=2:im
        fluxva(i,j)=fluxva(i,j)     ...
            -d(i,j)*2*aam2d(i,j)*(vab(i,j+1)-vab(i,j))     ...
            /dy(i,j);
    end
end
%
for j=2:jm
    for i=2:im
        fluxva(i,j)=fluxva(i,j)*dx(i,j);
        fluxua(i,j)=(fluxua(i,j)-tps(i,j))*0.25     ...
            *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1));
    end
end
%
for j=2:jmm1
    for i=2:imm1
        advva(i,j)=fluxua(i+1,j)-fluxua(i,j)     ...
            +fluxva(i,j)-fluxva(i,j-1);
    end
end
%
if(mode==2)
    %
    for j=2:jmm1
        for i=2:imm1
            wubot(i,j)=-0.50*(cbc(i,j)+cbc(i-1,j))     ...
                *sqrt(uab(i,j)^2     ...
                +(0.25*(vab(i,j)+vab(i,j+1)     ...
                +vab(i-1,j)+vab(i-1,j+1)))^2)     ...
                *uab(i,j);
        end
    end
    %
    for j=2:jmm1
        for i=2:imm1
            wvbot(i,j)=-0.50*(cbc(i,j)+cbc(i,j-1))     ...
                *sqrt(vab(i,j)^2     ...
                +(0.25*(uab(i,j)+uab(i+1,j)     ...
                +uab(i,j-1)+uab(i+1,j-1)))^2)     ...
                *vab(i,j);
        end
    end
    %
    for j=2:jmm1
        for i=2:imm1
            curv2d(i,j)=0.25     ...
                *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))     ...
                -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))     ...
                /(dx(i,j)*dy(i,j));
        end
    end
    %
    for j=2:jmm1
        for i=3:imm1
            advua(i,j)=advua(i,j)-aru(i,j)*0.25     ...
                *(curv2d(i,j)*d(i,j)     ...
                *(va(i,j+1)+va(i,j))     ...
                +curv2d(i-1,j)*d(i-1,j)     ...
                *(va(i-1,j+1)+va(i-1,j)));
        end
    end
    %
    for j=3:jmm1
        for i=2:imm1
            advva(i,j)=advva(i,j)+arv(i,j)*0.25     ...
                *(curv2d(i,j)*d(i,j)     ...
                *(ua(i+1,j)+ua(i,j))     ...
                +curv2d(i,j-1)*d(i,j-1)     ...
                *(ua(i+1,j-1)+ua(i,j-1)));
        end
    end
    %
end

return 


