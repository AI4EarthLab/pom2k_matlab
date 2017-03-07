function r = AXB(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    idx = 2:mx;
    A(idx,:,:) = 0.5*(A(idx-1,:,:) + A(idx,:,:));

    if(~isnumeric(obj_field))
%         switch obj_field.pos
%             case {0,2,4,6}
%                 A(1 ,: ,: )=0;  A(mx,: ,: )=0;
%             case {1,3,5,7}
%                 
%             otherwise
%                 disp('Unknown position.');
%         end
        
        A(1 ,: ,: )=0;
        r = Field(A, obj_field.grid, obj_field.grid.x_map(obj_field.pos)); 
    else
        r = A;
    end
end