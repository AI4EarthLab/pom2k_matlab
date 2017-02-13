function r = AZB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    idx = 2:kz;
    A(:,:,idx) = 0.5*(A(:,:,idx-1) + A(:,:,idx));
    
    if(~isnumeric(obj_field))
        switch obj_field.pos
            case {0,1,2,3}
                 
            case {4,5,6,7}
                A(: ,: ,1 )=0;  A(: ,: ,kz)=0;
            otherwise
                disp('Unknown position.');
        end 
%         % choose "constant" boundary condition or "zero gradient" boundary condition
%         switch obj_field.pos
%             case {0,1,2,3}
%                 A=A;
%             case {4,5,6,7}
%                 A(: ,: ,1 )=0;  A(: ,: ,kz)=0;
%             otherwise
%                 disp('Unknown position.');
%         end      
        r = Field(A, obj_field.grid, obj_field.grid.z_map(obj_field.pos));
    else
        r = A;
    end
end