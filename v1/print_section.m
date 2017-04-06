function print_section(iint,vamax,time,iext,tb,sb,dt_3d,et,imax,jmax)
global kbm1 dz_3d dx dy dx_3d dy_3d fsm_3d fsm iswtch prtd2 dti vmaxl sbias iprint

    if(iint>=iswtch)
        iprint=floor(prtd2*24.e0*3600.e0/dti + 0.5);
    end
    if(mod(iint,iprint)==0 || vamax>=vmaxl)
        fprintf('**************************************************\n');
        fprintf('**************************************************\n');
        fprintf('time = %.6f, iint = %d, iexit = %d, iprint = %d\n',time,iint,iext,iprint);
                
        tmp_dvol = dx_3d(:,:,1:kbm1).*dy_3d(:,:,1:kbm1).*fsm_3d(:,:,1:kbm1).*dt_3d(:,:,1:kbm1).*dz_3d(:,:,1:kbm1);
        vtot=sum(tmp_dvol(:));
        taver=sum(reshape(tb(:,:,1:kbm1).*tmp_dvol,[],1));
        saver=sum(reshape(sb(:,:,1:kbm1).*tmp_dvol,[],1));
        
        darea = dx.*dy.*fsm;    atot = sum(darea(:));
        eaver = sum(reshape(et.*darea,[],1));       
        taver=taver/vtot;       saver=saver/vtot;
        eaver=eaver/atot;       tsalt=(saver+sbias)*vtot;      
        fprintf('vtot = %.7f,atot = %.7f\n',vtot,atot);
        fprintf('eaver = %.7f,taver = %.7f,saver=%.7f,tsalt = %.7f\n',eaver,taver,saver,tsalt);
                
        if(vamax>vmaxl)            
            fprintf('time = %.6f, iint = %.6f, iexit = %.6f, iprint = %.6f\n',time,iint,iext,iprint);           
            disp        '************************************************'
            disp        '************ abnormal job end ******************'
            disp        '************* user terminated ******************'
            disp        '************************************************'
            
            fprintf('vamax = %d, imax = %d ,jmax = %d \n',vamax,imax,jmax);
            return;
        end
    end
return
