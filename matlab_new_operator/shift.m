function r = shift(A, unit, dim)
    b = circshift(A, unit, dim);

    if(dim == 1)
        if(unit > 0)
            b(1:unit, :, :) = A(1:unit, :, :);
        elseif(unit < 0)
            b(end+unit:end,:,:) = A(end+unit:end,:,:); 
        end
    end
    
    if(dim == 2)
        if(unit > 0)
            b(:, 1:unit, :) = A(:, 1:unit, :);
        elseif(unit < 0)
                b(:, end+unit:end,:) = A(:,end+unit:end,:); 
        end
    end
    
    if(dim == 3)
        if(unit > 0)
            b(:, :, 1:unit) = A(:, :, 1:unit);
        elseif(unit < 0)
            b(:, :, end+unit:end) = A(:, :, end+unit:end); 
        end
    end
    
    r = b;
end