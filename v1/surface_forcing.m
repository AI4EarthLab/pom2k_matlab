function [wusurf,wvsurf,wtsurf,wssurf,vfluxf,swrad,w] = surface_forcing(iint,time,iend,flag,t,s,vfluxf,wusurf,wvsurf,wtsurf,wssurf,swrad,w)
global im jm tbias sbias dti;
switch flag
    case {1,2,3}
        if(iint==1)
            vfluxf(:,:)=0.e0;   wtsurf(:,:)=0.0;
            satm=0.e0;          tatm=zeros(im,jm);
            tatm(:,:)=t(:,:,1)+tbias ;
            wtsurf=  wtsurf+vfluxf.*(tatm-t(:,:,1)-tbias);
            wssurf=  vfluxf.*(satm-s(:,:,1)-sbias)  ;
        else
            w(:,:,1)=vfluxf;
            return
        end
        
        %     case {2}
        %
        %     case {3}
        
    case {4}
    %============================surface:wind================================================    
        if wind_flag
            twind=30;                                     % time between wind files (days)
            iwind=int32(twind*86400.e0/dti);              % number of steps during a wind file
            % read wind stress data
            % read initial wind file
            if (iint == 1)
                [wusurff,wvsurff]=read_wind_pnetcdf(iint,iwind);
            end
            % read wind file corresponding to next time
            if (iint == 1 || mod(iint,iwind) == 0)
                wusurfb=wusurff;  wvsurfb=wvsurff;
                if (iint ~= iend)
                    [wusurff,wvsurff]=read_wind_pnetcdf(iint+iwind,iwind);
                end
            end
            % linear interpolation in time
            ntime=int32(time/twind);
            fnew=time/twind-ntime;
            fold=1.0-fnew;
            wusurf=fold*wusurfb+fnew*wusurff;
            wvsurf=fold*wvsurfb+fnew*wvsurff;
        end
    %=========================surface:heat====================================================    
        if heat_flag
            theat=1;                                %time between heat files(days)
            iheat=int32(theat*86400.e0/dti);        %number of steps during a heat file
            %read heat flux data
            % read initial heat file
            if (iint == 1)
                [wtsurff,swradf]=read_heat_pnetcdf(iint,iheat);
            end
            % read heat file corresponding to next time
            if (iint == 1 || mod(iint,iwind) == 0)
                wtsurfb=wtsurff;  swradb=swradf;
                if (iint ~= iend)
                    [wtsurff,swradf]=read_heat_pnetcdf(iint+iheat,iheat);
                end
            end
            % linear interpolation in time
            ntime=int32(time/theat);
            fnew=time/theat-ntime;
            fold=1.0-fnew;
            wtsurf=fold*wtsurfb+fnew*wtsurff;
            swrad=fold*swradb+fnew*swradf;
        end
    %====================surface:water===========================================================    
        if water_flag
            twater=1;                               %time between water files(days)
            iwater=int32(twater*86400.e0/dti);      %number of steps during a water file
            %read water flux data
            % read initial water file
            if (iint == 1)
                [wssurff]=read_water_pnetcdf(iint,iwater);
            end
            % read water file corresponding to next time
            if (iint == 1 || mod(iint,iwater) == 0)
                wssurfb=wssurff;
                if (iint ~= iend)
                    [wssurff]=read_heat_pnetcdf(iint+iwater,iwater);
                end
            end
            % linear interpolation in time
            ntime=int32(time/twater);
            fnew=time/twater-ntime;
            fold=1.0-fnew;
            wssurf=fold*wssurfb+fnew*wssurff;
        end
end
end