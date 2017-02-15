function obj = init_grid(obj, dx, dy, dz)

    obj.dx_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.dy_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.dz_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    % set point map relations        
    obj.x_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.y_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.z_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');    
    
    if(obj.type == 'C') 
        
        % set dx, dy, dz for each point, i.e. point 0, 1, 2, 3, 7
        obj.dx_map(0) = AYB(AXB(dx));
        obj.dy_map(0) = AYB(AXB(dy));
        obj.dz_map(0) = dz;
        
        obj.dx_map(1) = AYB(dx);
        obj.dy_map(1) = AYB(dy);
        obj.dz_map(1) = dz;
        
        obj.dx_map(2) = AXB(dx);
        obj.dy_map(2) = AXB(dy);
        obj.dz_map(2) = dz;

        obj.dx_map(3) = dx;
        obj.dy_map(3) = dy;
        obj.dz_map(3) = dz;
 
        obj.dx_map(4) = AYB(AXB(dx));
        obj.dy_map(4) = AYB(AXB(dy));
        obj.dz_map(4) = AZB(dz);

        obj.dx_map(5) = AYB(dx);
        obj.dy_map(5) = AYB(dy);
        obj.dz_map(5) = AZB(dz);

        obj.dx_map(6) = AXB(dx);
        obj.dy_map(6) = AXB(dy);
        obj.dz_map(6) = AZB(dz);
        
        obj.dx_map(7) = dx;
        obj.dy_map(7) = dy;
        obj.dz_map(7) = AZB(dz);

        obj.x_map(0) = 1;  obj.y_map(0) = 2;  obj.z_map(0) = 4;
        obj.x_map(1) = 0;  obj.y_map(1) = 3;  obj.z_map(1) = 5;
        obj.x_map(2) = 3;  obj.y_map(2) = 0;  obj.z_map(2) = 6;
        obj.x_map(3) = 2;  obj.y_map(3) = 1;  obj.z_map(3) = 7;
        obj.x_map(4) = 5;  obj.y_map(4) = 6;  obj.z_map(4) = 0;
        obj.x_map(5) = 4;  obj.y_map(5) = 7;  obj.z_map(5) = 1;    
        obj.x_map(6) = 7;  obj.y_map(6) = 4;  obj.z_map(6) = 2;
        obj.x_map(7) = 6;  obj.y_map(7) = 5;  obj.z_map(7) = 3;    
    end
    
    obj.has_initialized = true;
end