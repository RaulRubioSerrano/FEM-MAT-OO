classdef Triangle_Linear<Isoparametric
    %Triangle_Linear Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = Triangle_Linear
            obj = obj@Isoparametric;
            obj.type = 'TRIANGLE';
            obj.ndime = 2;
            obj.nnode = 3;
            obj.ngaus = 1;
            obj.weigp = 1/2;
            obj.posgp = [1/3;1/3];
            
            % s : xi coordinate
            % t : eta coordinate
            % u : zeta coordinate (for 3D)
            s = obj.posgp(1,obj.ngaus);
            t = obj.posgp(2,obj.ngaus);
            obj.deriv = zeros(obj.ndime,obj.nnode);
            obj.shape = zeros(1,obj.nnode);
            
            % Shape Functions
            obj.shape(1) = 1.0-s-t;
            obj.shape(2) = s;
            obj.shape(3) = t;
            
            %  SH Derivatives
            % w.r.t. xi
            obj.deriv(1,1) = -1.0;
            obj.deriv(1,2) = 1.0;
            obj.deriv(1,3) = 0.0;
            % w.r.t. eta
            obj.deriv(2,1) = -1.0;
            obj.deriv(2,2) = 0.0;
            obj.deriv(2,3) = 1.0;
        end    
    end
   
end