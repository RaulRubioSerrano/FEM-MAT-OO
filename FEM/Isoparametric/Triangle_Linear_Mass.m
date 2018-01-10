classdef Triangle_Linear_Mass<Isoparametric
    %Triangle_Linear Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = Triangle_Linear_Mass
            obj = obj@Isoparametric;
            obj.type = 'TRIANGLE';
            obj.ndime = 2;
            obj.nnode=3;
            obj.ngaus=3;
            obj.posgp(1,1)= 2.0/3.0;
            obj.posgp(2,1)= 1.0/6.0;
            obj.posgp(1,2)= 1.0/6.0;
            obj.posgp(2,2)= 2.0/3.0;
            obj.posgp(1,3)= 1.0/6.0;
            obj.posgp(2,3)= 1.0/6.0;
            obj.weigp(  1)= 1.0/6.0;
            obj.weigp(  2)= 1.0/6.0;
            obj.weigp(  3)= 1.0/6.0;
            
            obj.deriv = zeros(obj.ndime,obj.nnode,obj.ngaus);
            obj.shape = zeros(1,obj.nnode,obj.ngaus);
            for igauss=1:obj.ngaus
            % s : xi coordinate
            % t : eta coordinate
            % u : zeta coordinate (for 3D)
            s = obj.posgp(1,igauss);
            t = obj.posgp(2,igauss);

            % Shape Functions
            obj.shape(1,igauss) = 1.0-s-t;
            obj.shape(2,igauss) = s;
            obj.shape(3,igauss) = t;
            
            %  SH Derivatives
            % w.r.t. xi
            obj.deriv(1,1,igauss) = -1.0;
            obj.deriv(1,2,igauss) = 1.0;
            obj.deriv(1,3,igauss) = 0.0;
            % w.r.t. eta
            obj.deriv(2,1,igauss) = -1.0;
            obj.deriv(2,2,igauss) = 0.0;
            obj.deriv(2,3,igauss) = 1.0;
            end
        end    
    end
   
end