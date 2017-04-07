function [tbe,tbw,tbn,tbs,sbe,sbw,sbn,sbs,uabe,uabw,vabn,vabs] = lateral_bc(bc_flag,iint,time,iend,tbe,tbw,tbn,tbs,sbe,sbw,sbn,sbs,uabe,uabw,vabn,vabs)
global dti;
if(bc_flag)
      tbc=30; % time between bc files (days)
      ibc=int32(tbc*86400.e0/dti);
      ntime=time/tbc;
      % read bc data
      % read initial bc file
      if (iint==1)
       [tbwf,sbwf,uabwf,tbef,sbef,uabef,tbnf,sbnf,vabnf,tbsf,sbsf,vabsf]= ...
                 read_boundary_conditions_pnetcdf(iint,ibc);
      end
      % read bc file corresponding to next time
      if (iint == 1 || mod(iint,ibc) == 0)
            tbwb=tbwf;     sbwb=sbwf;
            tbeb=tbef;     sbeb=sbef;
            tbnb=tbnf;     sbnb=sbnf;
            tbsb=tbsf;     sbsb=sbsf;
          vabnb=vabnf;   vabsb=vabsf;
          uabwb=uabwf;   uabeb=uabef;
          if (iint ~= iend)
          [tbwf,sbwf,uabwf,tbef,sbef,uabef,tbnf,sbnf,vabnf,tbsf,sbsf,vabsf]= ...
                 read_boundary_conditions_pnetcdf(iint+ibc,ibc);        
          end
      end
% linear interpolation in time
      fnew=time/tbc-ntime;
      fold=1.e0-fnew;
      tbw=fold*tbwb+fnew*tbwf;
      sbw=fold*sbwb+fnew*sbwf;
      tbe=fold*tbeb+fnew*tbef;
      sbe=fold*sbeb+fnew*sbef;
      uabe=fold*uabeb+fnew*uabef;
      uabw=fold*uabwb+fnew*uabwf;
      tbn=fold*tbnb+fnew*tbnf;
      sbn=fold*sbnb+fnew*sbnf;
      tbs=fold*tbsb+fnew*tbsf;
      sbs=fold*sbsb+fnew*sbsf;
      vabn=fold*vabnb+fnew*vabnf;
      vabs=fold*vabsb+fnew*vabsf;
end
return     
end