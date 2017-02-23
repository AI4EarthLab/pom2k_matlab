function r = DZB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    A(:,:,2:kz) = A(:,:,2:kz) - A(:,:,1:kz-1);
    %A(:,:,1) = 0;
    if(~isnumeric(obj_field))
        % choose "constant" boundary condition or "zero gradient" boundary condition
        A(: ,: ,1 )=0;  A(: ,: ,kz)=0; 
        r = Field(A, obj_field.grid, obj_field.grid.z_map(obj_field.pos));
        dist=r.dz();
        r = DIVISION(r , dist(:,:,1:kz));         
    else
        r = A;
    end
end