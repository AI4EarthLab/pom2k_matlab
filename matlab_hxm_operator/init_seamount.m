function h=init_seamount(east_c,north_c)
    global im jm
    delh=0.9;
    %     Radius island or seamount:
    ra=25000.e0;
    %     Define depth:
    h=zeros(im,jm);
    for i=1:im
        for j=1:jm
            h(i,j)=4500.e0*(1.e0-delh *exp(-((east_c(i,j)...
                                        -east_c((im+1)/2,j))^2 ...
                                        +(north_c(i,j)...
                                         -north_c(i,(jm+1)/2))^2)...
                                         /ra^2));
            if(h(i,j) < 1.e0) 
                   h(i,j)=1.e0;
            end
        end 
    end 

    %     Close the north and south boundaries to form a channel:
    for i=1:im
        h(i,1)=1.e0;
        h(i,jm)=1.e0;
    end 
end