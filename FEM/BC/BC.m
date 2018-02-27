classdef BC<handle
    %BC Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties 
       
    end
    
    methods (Static,Access = public)
        % Constructor
        function bc = create(ptype) 
            switch ptype
                case {'ELASTIC','THERMAL'}
                    bc = BC_mechanics;
                case 'Stokes'
                    bc = BC_stokes;
            end
        end
    end
    
end

