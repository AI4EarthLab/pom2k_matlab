function r = DYB(obj_field)
    A = double(obj_field);
    dim = size(A);
    
    A(:,2:dim(2),:) = A(:,2:dim(2),:) - A(:,1:dim(2)-1,:);
    A(:,1,:) = 0;

    if(~isnumeric(obj_field))
        A = A ./ obj_field.dy_b();
        
        r = Field(A, obj_field.grid, ...
                  obj_field.grid.dyb_map(obj_field.pos));   
    else
        r = A;
    end
end