function r = DYB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    A(:,2:ny,:) = A(:,2:ny,:) - A(:,1:ny-1,:);

    if(~isnumeric(obj_field))
        % choose "constant" boundary condition or "zero gradient" boundary condition
        switch obj_field.pos
            case {0,1,4,5}              
                A(: ,1 ,: )=0;  A(: ,ny,: )=0; 
            case {2,3,6,7}
                 
            otherwise
                disp('Unknown position.');
        end    
        r = Field(A, obj_field.grid, obj_field.grid.y_map(obj_field.pos));  
        dist=r.dy();
        r = r ./ dist(:,:,1:kz);             
    else
        r = A;
    end
end