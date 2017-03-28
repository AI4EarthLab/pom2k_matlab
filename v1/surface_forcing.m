function [wusurf,wvsurf,wtsurf,wssurf,vfluxf,w] = surface_forcing(iint,time,iend,flag,t,s,vfluxf,wusurf,wvsurf,wtsurf,wssurf,w)
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

end