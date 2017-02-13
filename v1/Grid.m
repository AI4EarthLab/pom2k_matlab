classdef Grid
    properties
        dim
        type
        
        dx_map
        dy_map
        dz_map

        x_map
        y_map
        z_map
        
        has_initialized
    end
    methods
        function obj = Grid(dim_, type_)
            obj.dim = dim_;
            obj.type = type_;
            obj.has_initialized = false;
        end
        
        %     function obj = jump(g1, gt)
        %       t = mod(g1.type + gt, 2);
        %       obj = Grid(g1.dim, t, g1.gridset);
        %     end
        
        function r = dx(g1, pos)
            r = g1.dx_map(pos);
        end
        function r = dy(g1, pos)
            r = g1.dy_map(pos);
        end
        function r = dz(g1, pos)
            r = g1.dz_map(pos);
        end
        
    end
end

