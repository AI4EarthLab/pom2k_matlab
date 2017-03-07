    %
    %     Beginning of print section:
    %
    if(iint>=iswtch)
        iprint=floor(prtd2*24.e0*3600.e0/dti + 0.5);
    end
    %
    if(mod(iint,iprint)==0 || vamax>=vmaxl)
        fprintf('**************************************************\n');
        fprintf('**************************************************\n');
        fprintf('time = %.6f, iint = %d, iexit = %d, iprint = %d\n',time,iint,iext,iprint);
        
        %
        %     Select print statements in printall as desired:
        %
        %printall
        %
        vtot=0.e0;
        atot=0.e0;
        taver=0.e0;
        saver=0.e0;
        eaver=0.e0;
        
        dt_3d1 = repmat(dt, 1, 1, kb);
        tmp_dvol = dx_3d(:,:,1:kbm1).*dy_3d(:,:,1:kbm1).*fsm_3d(:,:,1:kbm1).*dt_3d1(:,:,1:kbm1).*dz_3d(:,:,1:kbm1);
        vtot=sum(tmp_dvol(:));
        taver=sum(reshape(tb(:,:,1:kbm1).*tmp_dvol,[],1));
        saver=sum(reshape(sb(:,:,1:kbm1).*tmp_dvol,[],1));
        
        %
%         for j=1:jm
%             for i=1:im
%                 darea=dx(i,j)*dy(i,j)*fsm(i,j);
%                 atot=atot+darea;
%                 eaver=eaver+et(i,j)*darea;
%             end
%         end
        %
        darea = dx.*dy.*fsm;
        atot = sum(darea(:));
        eaver = sum(reshape(et.*darea,[],1));
        
        taver=taver/vtot;
        saver=saver/vtot;
        eaver=eaver/atot;
        tsalt=(saver+sbias)*vtot;
      
        fprintf('vtot = %.6f,atot = %.6f\n',vtot,atot);
        fprintf('eaver = %.6f,taver = %.6f,saver=%.6f,tsalt = %.6f\n',eaver,taver,saver,tsalt);
        
        
        if(vamax>vmaxl)
            
            fprintf('time = %.6f, iint = %.6f, iexit = %.6f, iprint = %.6f\n',time,iint,iext,iprint);
           
%             printall(im,jm,imm1,jmm1,iskp,jskp,uab,vab,elb,d,dx,dy,time,u,v,w,t,s,rho,aam,km,kb,mode,dt,zz,z);
           
            disp        '************************************************'
            disp        '************ abnormal job end ******************'
            disp        '************* user terminated ******************'
            disp        '************************************************'
            
            fprintf('vamax = %d, imax = %d ,jmax = %d \n',vamax,imax,jmax);
            return;
            %
        end
        %
    end
    %
    %     End of print section