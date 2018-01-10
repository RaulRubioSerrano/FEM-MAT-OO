classdef DIM
    %DIM Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess = public, SetAccess = private)
        nstre
%         nnode
    end
    properties (GetAccess = {?Physical_Problem,?PhysicalVariables,?Element,?Postprocess}, SetAccess = private)
        ndim
%         ngaus
        nunkn
    end
    
    methods (Access = ?Physical_Problem)
        function obj = DIM(ptype,pdim)
            switch ptype
                case 'ELASTIC'
                    switch pdim
                        case '2D'
                            obj.ndim = 2;
                            obj.nunkn = 2;
                            obj.nstre = 3;
                        case '3D'
                            obj.ndim = 3;
                            obj.nunkn = 3;
                            obj.nstre = 6;
                    end
                case 'THERMAL'
                    error('Still not implemented.')
                case 'Stokes'
                        case '2D'
                            obj.ndim = 2;
                            obj.nunkn.u = 2;
                            obj.nunkn.p = 1;
                            obj.nstre = 0;
                        case '3D'
                            obj.ndim = 3;
                            obj.nunkn = 3;
                            obj.nstre = 6;
            end
        end
    end
    
end

