function gs=init_grid(gridtype)
global im jm kb dx_3d dy_3d dz_3d;
    gs = Grid([im,jm,kb], gridtype);
    gs.dx_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    gs.dy_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    gs.dz_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    % set point map relations        
    gs.x_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    gs.y_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    gs.z_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');    
    
    if(gs.type == 'C') 
        
        % set dx, dy, dz for each point, i.e. point 0, 1, 2, 3, 7
        gs.dx_map(0) = AYB(AXB(dx_3d));
        gs.dy_map(0) = AYB(AXB(dy_3d));
        gs.dz_map(0) = dz_3d;
        
        gs.dx_map(1) = AYB(dx_3d);
        gs.dy_map(1) = AYB(dy_3d);
        gs.dz_map(1) = dz_3d;
        
        gs.dx_map(2) = AXB(dx_3d);
        gs.dy_map(2) = AXB(dy_3d);
        gs.dz_map(2) = dz_3d;

        gs.dx_map(3) = dx_3d;
        gs.dy_map(3) = dy_3d;
        gs.dz_map(3) = dz_3d;
 
        gs.dx_map(4) = AYB(AXB(dx_3d));
        gs.dy_map(4) = AYB(AXB(dy_3d));
        gs.dz_map(4) = AZB(dz_3d);

        gs.dx_map(5) = AYB(dx_3d);
        gs.dy_map(5) = AYB(dy_3d);
        gs.dz_map(5) = AZB(dz_3d);

        gs.dx_map(6) = AXB(dx_3d);
        gs.dy_map(6) = AXB(dy_3d);
        gs.dz_map(6) = AZB(dz_3d);
        
        gs.dx_map(7) = dx_3d;
        gs.dy_map(7) = dy_3d;
        gs.dz_map(7) = AZB(dz_3d);

        gs.x_map(0) = 1;  gs.y_map(0) = 2;  gs.z_map(0) = 4;
        gs.x_map(1) = 0;  gs.y_map(1) = 3;  gs.z_map(1) = 5;
        gs.x_map(2) = 3;  gs.y_map(2) = 0;  gs.z_map(2) = 6;
        gs.x_map(3) = 2;  gs.y_map(3) = 1;  gs.z_map(3) = 7;
        gs.x_map(4) = 5;  gs.y_map(4) = 6;  gs.z_map(4) = 0;
        gs.x_map(5) = 4;  gs.y_map(5) = 7;  gs.z_map(5) = 1;    
        gs.x_map(6) = 7;  gs.y_map(6) = 4;  gs.z_map(6) = 2;
        gs.x_map(7) = 6;  gs.y_map(7) = 5;  gs.z_map(7) = 3;    
    end
    
    gs.has_initialized = true;
end