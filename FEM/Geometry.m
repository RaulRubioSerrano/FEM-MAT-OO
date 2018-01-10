classdef Geometry
    %Geometry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        cartDeriv
        shape
        dvolu
        weigp
        ndime
        nnode
        ngaus
        djacob
    end
    properties (GetAccess = ?Postprocess, SetAccess = private)
        posgp
    end
    
    methods (Access = ?Physical_Problem)
        function obj = Geometry(interpolation,quadrature,nelem)
%             obj.nnode = size(mesh.connec,2);
%             switch mesh.geometryType
%                 case 'TRIANGLE'
%                     switch obj.nnode
%                         case 3
%                             geometryObject = Triangle_Linear;
%                         case 6
%                             geometryObject = Triangle_Quadratic;
%                         otherwise
%                             error('Invalid nnode for element TRIANGLE.');
%                     end
%                 case 'Triangle_Linear_Mass'
%                     geometryObject=Triangle_Linear_Mass;
%                 case 'QUAD'
%                     switch obj.nnode
%                         case 4
%                             geometryObject = Quadrilateral_Bilinear;
%                         case 8
%                             geometryObject = Quadrilateral_Serendipity;
%                         otherwise
%                             error('Invalid nnode for element QUADRILATERAL.');
%                     end
%                 case 'TETRAHEDRA'
%                     geometryObject = Tetrahedra;
%                 case 'HEXAHEDRA'
%                     geometryObject = Hexahedra;
%                 otherwise
%                     error('Invalid mesh type.')
%             end
%           
            
            geometryObject = interpolation.isoparametric;
            obj.ndime = interpolation.isoparametric.ndime;
            obj.posgp = quadrature.posgp;
            obj.weigp = quadrature.weigp;
            obj.ngaus = quadrature.ngaus;
%             obj.shape = interpolation.isoparametric.shape;


            
            
            for i = 1:obj.ndime
                a = interpolation.xpoints(:,i);
                elcoord(:,i,:) = a(permute(interpolation.T',[1,3,2]));
            end
            % Gauss loop
            for igauss = 1:obj.ngaus
                 pos_gp = num2cell(obj.posgp(1:obj.ndime,igauss));
                 obj.shape(:,igauss) = cell2mat(interpolation.isoparametric.shape(pos_gp{:}));
                for i = 1:obj.ndime
                    % Jacobian
                    deriv_perm = permute(geometryObject.deriv(i,:,igauss),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,obj.ndime,nelem);
                    jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
                end                
                % !! SWITCH EXECPCIÓ? !!
                switch obj.ndime
                    case 2
                        [invJ,detJ] = multinverse2x2(jacobian);
                    case 3
                        [invJ,detJ] = multinverse3x3(jacobian);
                end
                
                for i = 1:obj.ndime
                    % Cartesian Derivatives
                    deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,geometryObject.nnode,1) .*repmat(geometryObject.deriv(:,:,igauss),1,1,nelem);
                    obj.cartDeriv(i,:,:,igauss) = sum(deriv_perm_large);
                end
                obj.dvolu(:,igauss) = geometryObject.weigp(igauss)*detJ;
                obj.djacob(:,igauss)= detJ;
            end
            % dvolu
            
            
        end
    end
    
end

