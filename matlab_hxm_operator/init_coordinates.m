function [dx,dy,dz,z,zz,dzz,east_c,north_c,rot]=init_coordinates()

    [z,zz,dz,dzz]       = init_depth();
    [dx,dy]             = init_latlon();
    [east_c,north_c,rot]= init_horizon_rot(dx,dy);

    
end


function [z,zz,dz,dzz]=init_depth()
    global kb kl1 kl2
    kdz=[1,1,2,4,8,16,32,64,128,256,512,1024];

    z=zeros(1,kb);
    zz=zeros(1,kb);
    dz=zeros(1,kb);
    dzz=zeros(1,kb);

    for k=2:kl1
        z(k)=z(k-1)+kdz(k-1);
    end

    delz=z(kl1)-z(kl1-1);

    for k=kl1+1:kl2
        z(k)=z(k-1)+delz;
    end

    for k=kl2+1:kb
        dz(k)=kdz(kb-k+1)*delz/kdz(kb-kl2);
        z(k)=z(k-1)+dz(k);
    end

    for k=1:kb
        z(k)=-z(k)/z(kb);
    end

    for k=1:kb-1
        zz(k)=0.5e0*(z(k)+z(k+1));
    end

    zz(kb)=2.e0*zz(kb-1)-zz(kb-2);

    for k=1:kb-1
        dz(k)=z(k)-z(k+1);
        dzz(k)=zz(k)-zz(k+1);
    end

    dz(kb)=0.e0;
    dzz(kb)=0.e0;
end 


function [dx,dy]=init_latlon()
    global im jm 
    global pi
    delx=8000.e0;

    for j=1:jm
        for i=1:im
        %     For constant grid size:
        %         dx(i,j)=delx
        %         dy(i,j)=delx
        %     Set up grid dimensions, areas of free surface cells, and
        %     Coriolis parameter:
        %
        %     For variable grid size:
            dx(i,j)=delx-delx*sin(pi*i/im)/2.0;
            dy(i,j)=delx-delx*sin(pi*j/jm)/2.0;
        end
    end
end


function [east_c,north_c,rot]=init_horizon_rot(dx,dy)
    global im jm 
    %     Calculate horizontal coordinates of grid points and rotation
    %     angle.
    %
    %     NOTE that this is introduced solely for the benefit of any post-
    %     processing software, and in order to conform with the requirements
    %     of the NetCDF Climate and Forecast (CF) Metadata Conventions.
    %
    %     There are four horizontal coordinate systems, denoted by the
    %     subscripts u, v, e and c ("u" is a u-point, "v" is a v-point,
    %     "e" is an elevation point and "c" is a cell corner), as shown
    %     below. In addition, "east_*" is an easting and "north_*" is a
    %     northing. Hence the coordinates of the "u" points are given by
    %     (east_u,north_u).
    %
    %     Also, if the centre point of the cell shown below is at
    %     (east_e(i,j),north_e(i,j)), then (east_u(i,j),north_u(i,j)) are
    %     the coordinates of the western of the two "u" points,
    %     (east_v(i,j),north_v(i,j)) are the coordinates of the southern of
    %     the two "v" points, and (east_c(i,j),north_c(i,j)) are the
    %     coordinates of the southwestern corner point of the cell. The
    %     southwest corner of the entire grid is at
    %     (east_c(1,1),north_c(1,1)).
    %
    %                      |              |
    %                    --c------v-------c--
    %                      |              |
    %                      |              |
    %                      |              |
    %                      |              |
    %                      u      e       u
    %                      |              |
    %                      |              |
    %                      |              |
    %                      |              |
    %                    --c------v-------c--
    %                      |              |
    %
    %
    %     NOTE that the following calculation of east_c and north_c only
    %     works properly for a rectangular grid with east and north aligned
    %     with i and j, respectively:
    %
    %    Compute east_c and north_c
    for j=1:jm
        east_c(1,j)=0.0;
        for i=2:im
            east_c(i,j)=east_c(i-1,j)+dx(i-1,j);
        end
    end

    for i=1:im
        north_c(i,1)=0.0;
        for j=2:jm
            north_c(i,j)=north_c(i,j-1)+dy(i,j-1);
        end
    end

    %
    %     The following works properly for any grid:
    %
    %     Elevation points:
    %
    for j=1:jm-1
        for i=1:im-1
            east_e(i,j)=(east_c(i,j)+east_c(i+1,j)+east_c(i,j+1)+east_c(i+1,j+1))/4.e0;
            north_e(i,j)=(north_c(i,j)+north_c(i+1,j)+north_c(i,j+1)+north_c(i+1,j+1))/4.e0;
        end
    end
    %     Extrapolate ends:
    for i=1:im-1
        east_e(i,jm)=((east_c(i,jm)+east_c(i+1,jm))*3.e0-east_c(i,jm-1)-east_c(i+1,jm-1))/4.e0;
        north_e(i,jm)=((north_c(i,jm)+north_c(i+1,jm))*3.e0-north_c(i,jm-1)-north_c(i+1,jm-1))/4.e0;
    end
    %
    for j=1:jm-1
        east_e(im,j)=((east_c(im,j)+east_c(im,j+1))*3.e0-east_c(im-1,j)-east_c(im-1,j+1))/4.e0;
        north_e(im,j)=((north_c(im,j)+north_c(im,j+1))*3.e0-north_c(im-1,j)-north_c(im-1,j+1))/4.e0;
    end
    %
    east_e(im,jm)=east_e(im-1,jm)+east_e(im,jm-1)-(east_e(im-2,jm)+east_e(im,jm-2))/2.e0;
    north_e(im,jm)=north_e(im-1,jm)+north_e(im,jm-1)-(north_e(im-2,jm)+north_e(im,jm-2))/2.e0;

    %     u-points:
    for j=1:jm-1
        for i=1:im
            east_u(i,j)=(east_c(i,j)+east_c(i,j+1))/2.e0;
            north_u(i,j)=(north_c(i,j)+north_c(i,j+1))/2.e0;
        end
    end
    %     Extrapolate ends:
    for i=1:im
        east_u(i,jm)=(east_c(i,jm)*3.0-east_c(i,jm-1))/2.e0;
        north_u(i,jm)=(north_c(i,jm)*3.0-north_c(i,jm-1))/2.e0;
    end

    %     v-points:
    for j=1:jm
        for i=1:im-1
            east_v(i,j)=(east_c(i,j)+east_c(i+1,j))/2.e0;
            north_v(i,j)=(north_c(i,j)+north_c(i+1,j))/2.e0;
        end
    end
    %     Extrapolate ends:
    for j=1:jm
        east_v(im,j)=(east_c(im,j)*3.0-east_c(im-1,j))/2.e0;
        north_v(im,j)=(north_c(im,j)*3.0-north_c(im-1,j))/2.e0;
    end

    %     rot is the angle (radians, anticlockwise) of the i-axis relative
    %     to east, averaged to a cell centre:
    %
    %     (NOTE that the following calculation of rot only works properly
    %     for this particular rectangular grid)
    %
    rot=zeros(im,jm);

end



