function r = shift(A, unit, dim)
    b = circshift(A.data, unit, dim);
    B=double(A);
    [mx,ny,kz] = size(B);
    if(dim == 1)
        if(unit > 0)
            b(1:unit, :, :) = A(1:unit, :, :);
        elseif(unit < 0)
            b(mx+unit:mx,:,:) = A(mx+unit:mx,:,:); 
        end
    end
    
    if(dim == 2)
        if(unit > 0)
            b(:, 1:unit, :) = A(:, 1:unit, :);
        elseif(unit < 0)
                b(:, ny+unit:ny,:) = A(:,ny+unit:ny,:); 
        end
    end
    
    if(dim == 3)
        if(unit > 0)
            b(:, :, 1:unit) = A(:, :, 1:unit);
        elseif(unit < 0)
            b(:, :, kz+unit:kz) = A(:, :, kz+unit:kz); 
        end
    end
    
    r = b;
end