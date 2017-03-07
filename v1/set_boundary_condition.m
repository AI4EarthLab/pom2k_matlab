
    %-----------------------------------------------------------------------
    %     Set time dependent, surface and lateral boundary conditions.
    %     The latter will be used in subroutine bcond. Users may
    %     wish to create a subroutine to supply wusurf, wvsurf, wtsurf,
    %     wssurf, swrad and vflux.
    %
    %     Introduce simple wind stress. Value is negative for westerly or
    %     southerly winds. The following wind stress has been tapered
    %     along the boundary to suppress numerically induced oscilations
    %     near the boundary (Jamart and Ozer, J.G.R., 91, 10621-10631).
    %     To make a healthy surface Ekman layer, it would be well to set
    %     kl1=9.
    %

for j=2:jmm1
        for i=2:imm1
            if(iproblem~=3) % constant wind read in file2ic
    %           wusurf(i,j)=ramp*(1.e-4*cos(pi*(j-1)/jmm1));
                wusurf(i,j)=1.00*(1.e-4*cos(pi*(j-1)/jmm1))  ...
                    *0.25*(dvm(i,j+1)+dvm(i-1,j+1)     ...
                    +dvm(i-1,j)+dvm(i,j));
                % --- no wind ----
                %           wusurf(i,j)=0.e0;
                wvsurf(i,j)=0.e0;
            end
            
            e_atmos(i,j)=0.e0;
            vfluxf(i,j)=0.e0;
            %
            %     Set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
            %     the sea surface. See calculation of elf(i,j) below and subroutines
            %     vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
            %     is no net flow across lateral boundaries, the basin volume will be
            %     constant; if also vflux(i,j).ne.0, then, for example, the average
            %     salinity will change and, unrealistically, so will total salt.
            %
            w(i,j,1)=vfluxf(i,j);
            %
            %     Set wtsurf to the sensible heat, the latent heat (which involves
            %     only the evaporative component of vflux) and the long wave
            %     radiation:
            %
            wtsurf(i,j)=0.0;
            %
            %     Set swrad to the short wave radiation:
            %
            swrad(i,j)=0.0;
            %
            %     To account for change in temperature of flow crossing the sea
            %     surface (generally quite small compared to latent heat effect)
            %
            tatm=t(i,j,1)+tbias;    % an approximation
            wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias);
            %
            %     Set the salinity of water vapor/precipitation which enters/leaves
            %     the atmosphere (or e.g., an ice cover)
            %
            satm=0.0  ;
            wssurf(i,j)=  vfluxf(i,j)*(satm-s(i,j,1)-sbias)  ;
            %
        end
    end