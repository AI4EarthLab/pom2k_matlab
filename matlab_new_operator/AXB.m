function r = AXB(obj_field)
    A = double(obj_field);
    dim = size(A);
    
    idx = 2:dim(1);
    A(idx,:,:) = 0.5*(A(idx-1,:,:) + A(idx,:,:));

    if(~isnumeric(obj_field))
        r = Field(A, obj_field.grid, ...
                  obj_field.grid.axb_map(obj_field.pos));
    else
        r = A;
    end
end