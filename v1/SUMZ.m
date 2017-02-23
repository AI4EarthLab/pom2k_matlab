function r = SUMZ(obj_field)
    A = double(obj_field);
    [mx,ny,kz] = size(A);
    
    B=zeros(mx,ny,kz);
    B(:,:,1)=A(:,:,1);
    for k=2:kz-1
        B(:,:,k)=B(:,:,k-1)+A(:,:,k);
    end   
    if(~isnumeric(obj_field))
        r = Field(B, obj_field.grid, obj_field.pos);     
    else
        r = B;
    end
end