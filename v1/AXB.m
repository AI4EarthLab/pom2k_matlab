function r = AXB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    idx = 2:mx;
    A(idx,:,:) = 0.5*(A(idx-1,:,:) + A(idx,:,:));

    if(~isnumeric(obj_field))        
        A(1 ,: ,: )=0;
        r = Field(A, obj_field.grid, obj_field.grid.x_map(obj_field.pos)); 
    else
        r = A;
    end
end