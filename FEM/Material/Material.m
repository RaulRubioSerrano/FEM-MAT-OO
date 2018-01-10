classdef Material
    %Material Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = protected, SetAccess = private)
        nelem
    end
    
    methods (Access = protected)
        function obj = Material(nelem)
            obj.nelem = nelem;
        end
    end
    
    methods (Access = ?Physical_Problem, Static)
        function material = create(ptype,pdim,nelem)
            switch ptype
                case 'ELASTIC'
                    % !! IT HAS BEEN ASSUMED THAT THERE'S ONLY ISOTROPIC MATERIALS. 
                    %    THIS HAS TO BE CHANGED FOR THE OPT TOP PROBLEM !!
                    switch pdim
                        case '2D'
                            material = Material_Elastic_ISO_2D(nelem);
                        case '3D'
                            material = Material_Elastic_ISO_3D(nelem);
                    end
                    
                case 'THERMAL'
                    error('Still not implemented.')
                case 'Stokes'
                    material = Material_Stokes(nelem);
                otherwise
                    error('Invalid ptype.')
            end
        end
    end
    
end

