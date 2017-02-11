function r = DXF(obj_field)
    A = double(obj_field);
    dim = size(A);
    
    A(1:dim(1)-1,:,:) = A(2:dim(1),:,:) - A(1:dim(1)-1,:,:);
    A(dim(1),:,:) = 0;

    if(~isnumeric(obj_field))
        A = A ./ obj_field.dx_f();
        
        r = Field(A, obj_field.grid, ...
                  obj_field.grid.dxf_map(obj_field.pos));   
    else
        r = A;
    end
end