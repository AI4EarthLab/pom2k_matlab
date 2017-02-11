function r = DZF(obj_field)
    A = double(obj_field);
    dim = size(A);
    
    idx = 1:dim(3)-1;
    
    A(:,:,idx) = A(:,:,idx+1) - A(:,:,idx);
    A(:,:,dim(3)) = 0;

    if(~isnumeric(obj_field))
        A = A ./ obj_field.dz_f();
        
        r = Field(A, obj_field.grid, ...
                  obj_field.grid.dzf_map(obj_field.pos));   
    else
        r = A;
    end
end