function r = AZF(obj_field)
    A = double(obj_field);
    dim = size(A);
    
    idx = 1:dim(3)-1;
    A(:,:,idx) = 0.5*(A(:,:,idx) + A(:,:,idx+1));
    
    if(~isnumeric(obj_field))
        r = Field(A, obj_field.grid, ...
                  obj_field.grid.azf_map(obj_field.pos));
    else
        r = A;
    end
end