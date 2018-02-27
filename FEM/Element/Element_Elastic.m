classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
%     methods (Access = ?Physical_Problem)
%         
%     end
    methods (Access = protected)
        function Ke = computeMatrix(obj,dim,nelem,geometry_test,geometry_variable,material,ifield,jfield)
           
           
            Ke = zeros(dim.nunkn*geometry_test.nnode,dim.nunkn*geometry_variable.nnode,nelem);
            % Elastic matrix
             Cmat = material.C;           
            
            for igauss = 1 :geometry_variable.ngaus
                % Strain-displacement matrix
                [obj.B, Bmat] = obj.B.computeB(dim.nunkn,nelem,geometry_variable.nnode,geometry_variable.cartDeriv(:,:,:,igauss));
                
                for iv=1:geometry_test.nnode*dim.nunkn
                    for jv=1:geometry_variable.nnode*dim.nunkn
                        for istre=1:dim.nstre
                            for jstre=1:dim.nstre
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*geometry_variable.dvolu(:,igauss);
                            end
                        end
                        
                    end
                end
            end
        end
        
        function Fext = computeSuperficialRHS(obj,dim,nelem,geometry_variable,bc,dof) %To be donne
            Fext = zeros(geometry_variable.nnode*dim.nunkn,1,nelem);
        end
        function Fext = computeVolumetricRHS(obj,dim,nelem,geometry_variable,bc,dof)%To be done
            Fext = zeros(geometry_variable.nnode*dim.nunkn,1,nelem);
            
        end
    end
    
end
