classdef Grid
    properties
        dim
        type
        
        dx_f_map
        dy_f_map
        dz_f_map

        dx_b_map
        dy_b_map
        dz_b_map
        
        axf_map
        axb_map
        ayf_map
        ayb_map
        azf_map
        azb_map
        
        dxf_map
        dxb_map
        dyf_map
        dyb_map
        dzf_map
        dzb_map

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
        
        function r = dx_f(g1, pos)
            r = g1.dx_f_map(pos);
        end
        function r = dy_f(g1, pos)
            r = g1.dy_f_map(pos);
        end
        function r = dz_(g1, pos)
            r = g1.dz_f_map(pos);
        end
        
        function r = dx_b(g1, pos)
            r = g1.dx_b_map(pos);
        end
        function r = dy_b(g1, pos)
            r = g1.dy_b_map(pos);
        end
        function r = dz_b(g1, pos)
            r = g1.dz_b_map(pos);
        end
        
    end
end

