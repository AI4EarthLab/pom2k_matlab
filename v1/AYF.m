function r = AYF(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    idx = 1:ny-1;
    A(:,idx,:) = 0.5*(A(:,idx,:) + A(:,idx+1,:));
    
    if(~isnumeric(obj_field))
        switch obj_field.pos
            case {0,1,4,5}
                A(: ,1 ,: )=0;  A(: ,ny,: )=0;
            case {2,3,6,7}
                
            otherwise
                disp('Unknown position.');
        end     
        r = Field(A, obj_field.grid, obj_field.grid.y_map(obj_field.pos)); 
    else
        r = A;
    end
end