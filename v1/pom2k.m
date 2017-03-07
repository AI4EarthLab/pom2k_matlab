clear all;
init_global();
init_constants();
init_operator();
create_grid();
init_variables();
read_initvalues();
init_grid();
init_fields();
init_time_condition();

%==========================================
%           begin internal (3-D) mode
%==========================================
for iint=1:iend
    time=dti*iint*1.0/86400.0+time0;
    if(lramp~=0)
        ramp = time/period;
        if(ramp>1.0)
            ramp=1.0;
        end
    else
        ramp=1.0;
    end

    set_boundary_condition();

    cal_adxy2d();  
    
    cal_external_mode();
    
    if(vamax<=vmaxl) 

        %     continue with internal (3-D) mode calculation:
        if((iint~= 1|| time0~=0.e0) && mode~=2)

            adjust_uv();
          
            w     = vertvl (dt_3d,u,v,vfluxb,vfluxf,etf,etb);  

            cal_q2q2l();       

            cal_ts();

            cal_uv();


        end  % end if

        update_3dvar();

    end   
    
    print_section();
    %-----------------------------------------------------------------------
    %
    %  End of internal (3-D) mode
    %-----------------------------------------------------------------------  
end

fprintf('time = %d, iint = %d ,iext = %d , iprint = %d \n',time,iint,iext,iprint);

