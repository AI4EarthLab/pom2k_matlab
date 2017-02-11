function r = DZB(obj_field)
    A = double(obj_field);
    dim = size(A)
    
    idx = 2:dim(3);
    
    A(:,:,idx) = A(:,:,idx) - A(:,:,idx-1);
    A(:,:,1) = 0;

    if(~isnumeric(obj_field))
        A = A ./ obj_field.dz_b();
        
        r = Field(A, obj_field.grid, ...
                  obj_field.grid.dzb_map(obj_field.pos));   
    else 
        r = A;
    end
end