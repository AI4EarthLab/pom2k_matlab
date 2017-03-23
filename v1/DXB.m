function r = DXB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    A(2:mx,:,:) = A(2:mx,:,:) - A(1:mx-1,:,:);

    if(~isnumeric(obj_field))
        % choose "constant" boundary condition or "zero gradient" boundary condition
        switch obj_field.pos
            case {0,2,4,6}
                A(1 ,: ,: )=0;  A(mx,: ,: )=0;
%             case {1,3,5,7}
%                  
%             otherwise
%                 disp('Unknown position.');
        end
        r = Field(A, obj_field.grid, obj_field.grid.x_map(obj_field.pos));   
        dist=r.dx();
        r = r ./ dist(:,:,1:kz);  
    else
        r = A;
    end
end