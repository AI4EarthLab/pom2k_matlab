function obj = init_grid(obj, dx, dy, dz)

    obj.dx_f_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.dy_f_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.dz_f_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');

    obj.dx_b_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.dy_b_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    obj.dz_b_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');

    if(obj.type == 'C') 
        
        % set dx, dy, dz for each point, i.e. point 0, 1, 2, 3, 7
        obj.dx_f_map(0) = AYB(AXB(dx));
        obj.dy_f_map(0) = AYB(AXB(dy));
        obj.dz_f_map(0) = AZF(dz);
        
        obj.dx_f_map(1) = AYB(dx);
        obj.dy_f_map(1) = AYB(dy);
        obj.dz_f_map(1) = AZF(dz);
        
        obj.dx_f_map(2) = AXB(dx);
        obj.dy_f_map(2) = AXB(dy);
        obj.dz_f_map(2) = AZF(dz);

        obj.dx_f_map(3) = dx;
        obj.dy_f_map(3) = dy;
        obj.dz_f_map(3) = AZF(dz);
        
        obj.dx_f_map(7) = dx;
        obj.dy_f_map(7) = dy;
        obj.dz_f_map(7) = dz;

        obj.dx_b_map(0) = shift(AYB(AXB(dx)), -1, 1);
        obj.dy_b_map(0) = shift(AYB(AXB(dy)), -1, 2);
        obj.dz_b_map(0) = shift(AZF(dz), -1, 3);
        
        obj.dx_b_map(1) = shift(AYB(dx), -1, 1);
        obj.dy_b_map(1) = shift(AYB(dy), -1, 2);
        obj.dz_b_map(1) = shift(AZF(dz), -1, 3);
        
        obj.dx_b_map(2) = shift(AXB(dx), -1, 1);
        obj.dy_b_map(2) = shift(AXB(dy), -1, 2);
        obj.dz_b_map(2) = shift(AZF(dz), -1, 3);

        obj.dx_b_map(3) = shift(dx,      -1, 1);
        obj.dy_b_map(3) = shift(dy,      -1, 2);
        obj.dz_b_map(3) = shift(AZF(dz), -1, 3);
        
        obj.dx_b_map(7) = shift(dx, -1, 1);
        obj.dy_b_map(7) = shift(dy, -1, 2);
        obj.dz_b_map(7) = shift(dz, -1, 3);
        
        % set point map relations        
        dxf_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        dxb_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        dyf_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        dyb_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        dzf_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');        
        dzb_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');

        axf_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        axb_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        ayf_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        ayb_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
        azf_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');        
        azb_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');

        dxf_map(0) = 1;  dyf_map(0) = 2;  dzf_map(0) = 4;
        dxf_map(1) = 0;  dyf_map(1) = 3;  dzf_map(1) = 5;
        dxf_map(2) = 3;  dyf_map(2) = 0;  dzf_map(2) = 6;
        dxf_map(3) = 2;  dyf_map(3) = 1;  dzf_map(3) = 7;
        dxf_map(4) = 5;  dyf_map(4) = 6;  dzf_map(4) = 0;
        dxf_map(5) = 4;  dyf_map(5) = 7;  dzf_map(5) = 1;    
        dxf_map(6) = 7;  dyf_map(6) = 4;  dzf_map(6) = 2;
        dxf_map(7) = 6;  dyf_map(7) = 5;  dzf_map(7) = 3;    

        axf_map = dxf_map;  ayf_map = dyf_map;  azf_map = dzf_map;
        axb_map = dxf_map;  ayb_map = dyf_map;  azb_map = dzf_map;
        dxb_map = dxf_map;  dyb_map = dyf_map;  dzb_map = dzf_map;
        
        obj.axf_map = axf_map;  obj.dxf_map = dxf_map;
        obj.axb_map = axb_map;  obj.dxb_map = dxb_map;
        obj.ayf_map = ayf_map;  obj.dyf_map = dyf_map;
        obj.ayb_map = ayb_map;  obj.dyb_map = dyb_map;
        obj.azf_map = azf_map;  obj.dzf_map = dzf_map;
        obj.azb_map = azb_map;  obj.dzb_map = dzb_map;
    end
    
    obj.has_initialized = true;
end