classdef Element_Thermal < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = {?Physical_Problem, ?Element})
        function obj = Element_Thermal()
            obj.nincr = 1;
            obj.cload = 0;
        end
        

        function [r,dr] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************
            [K] = obj.computeStiffnessMatrix();
                     
            R = obj.compute_imposed_displacemet_force(K);
            fext = obj.cload + R;
            
            fint = K(obj.dof.vF,obj.dof.vF)*uL;

            r = fint - fext;
            dr = K(obj.dof.vF,obj.dof.vF);
        end
        
        
        function K = computeStiffnessMatrix(obj)
            K = compute_elem_StiffnessMatrix(obj);
            K = obj.AssembleMatrix(K);
        end
        
        function K = compute_elem_StiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            
            
            for igaus = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                Bmat = obj.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                % Compute Ke
                if obj.nelem < 1000 %Just to reduce test.m compute time TO BE REMOVED
                    for i = 1:obj.nelem
                        Ke(:,:,i) = Ke(:,:,i)+Bmat(:,:,i)'*Bmat(:,:,i)*obj.geometry.dvolu(i,igaus);
                    end
                else
                    for iv = 1:obj.nnode*obj.nunkn
                        for jv = 1:obj.nnode*obj.nunkn
                            for istre = 1:obj.nstre
                                % for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igaus);
                                %end
                            end
                        end
                    end
                end
            end
            K = Ke;
        end
        
        function [B] = computeB(obj,nunkn,nelem,nnode,cartd)
            B = zeros(2,nnode*nunkn,nelem);
            for inode=1:nnode
                j = nunkn*(inode-1)+1;
                B(1,j,:)=cartd(1,inode,:);
                B(2,j,:)=cartd(2,inode,:);
            end
            
            
        end
    
    end

    
    methods (Access = protected)
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
    end
    
end


