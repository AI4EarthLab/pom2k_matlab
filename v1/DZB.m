function r = DZB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    A(:,:,2:kz) = A(:,:,2:kz) - A(:,:,1:kz-1);
    %A(:,:,1) = 0;
    if(~isnumeric(obj_field))
        % choose "constant" boundary condition or "zero gradient" boundary condition
        switch obj_field.pos
            case {0,1,2,3}
                 
            case {4,5,6,7}
                A(: ,: ,1 )=0;  A(: ,: ,kz)=0;
            otherwise
                disp('Unknown position.');
        end      
        r = Field(A, obj_field.grid, obj_field.grid.z_map(obj_field.pos));
        dist=r.dz();
        r = DIVISION(r , dist(:,:,1:kz));         
    else
        r = A;
    end
end