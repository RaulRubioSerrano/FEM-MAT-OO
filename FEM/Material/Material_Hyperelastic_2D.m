classdef Material_Hyperelastic_2D < Material_Elastic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Element_Hyperelastic,?PhysicalVars_Elastic, ?Physical_Problem}, SetAccess = ?Physical_Problem)
        connec
        cartd
        nnode
        coord
    end
    
    methods
        function obj = Material_Hyperelastic_2D(nelem,connec,cartd,nnode,coord)
            obj@Material_Elastic(nelem);
            obj.connec= connec;
            obj.cartd = cartd;
            obj.nnode = nnode;
            obj.coord = coord;
        end
        
        %% Compute Eulerian elasticity tensor
        function [ctens,sigma] = computeCtens(obj,coord)
            
            % Jacobian
            [F,Fjacb,blcg,Crcg] = obj.computeDefGradient(coord,obj.cartd);
            
            % Effective Lame moduli
            mup     = (obj.mu-obj.lambda*log(Fjacb))./Fjacb;
            lambdap = obj.lambda./Fjacb;
            
            % Rearrange the dimensions
            mup     = permute(mup,[2 3 1]);
            lambdap = permute(lambdap,[2 3 1]);
            
            % Cauchy stress tensor for a compressible neo-Hookean material
            sigma = obj.computeCauchyStress(F,Fjacb,blcg,Crcg); 
            
            % Define constitutive tensor
            ctens = zeros(3,3,3,3,obj.nelem); % material
            
            % Define delta kronecker
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem]);
            
            % Compute c
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            inc = lambdap.*dk(i,j,:).*dk(k,l,:) + mup.*(dk(i,k,:).*dk(j,l,:) + dk(i,l,:).*dk(j,k,:));
                            ctens(i,j,k,l,:) = squeeze(ctens(i,j,k,l,:)) + squeeze(inc);
                        end
                    end
                end
            end
            
        end
        
        %% Compute Cauchy stress tensor
        function sigma = computeCauchyStress(obj,F,Fjacb,blcg,Crcg)
            I = repmat(eye(3),[1 1 obj.nelem]);
            Fjacb = permute(Fjacb,[2 3 1]);

            mu      = repmat(obj.mu,[1 1 obj.nelem]);
            lambda  = repmat(obj.lambda,[1 1 obj.nelem]);

            % Define delta kronecker
            dk = eye(3); dk = repmat(dk,[1 1 obj.nelem]);
            
            % Strain
            E = 1/2*(Crcg-I);
            
            % Deformation tensor inverse
            Cinv = multinverse3x3(Crcg);
            
            % Second Piola-Kirchhoff stress tensor
            S = zeros(3,3,obj.nelem);
            
            for i = 1:3
                for j = 1:3
                    S(i,j,:) = squeeze(mu).*squeeze(I(i,j,:)-Cinv(i,j,:)) + squeeze(lambda).*squeeze(log(Fjacb)).*squeeze(Cinv(i,j,:));
                end
            end
            
%            St. Venant-Kirchhoff
%             for i = 1:3
%                 for j = 1:3
%                     for k = 1:3
%                     S(i,j,:) = S(i,j,:) + obj.lambda*E(k,k,:)*dk(i,j);
%                     end
%                     S(i,j,:) = S(i,j,:) + 2*obj.mu*E(i,j,:);
%                 end
%             end

            
            sigma = zeros(3,3,obj.nelem);
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            sigma(i,l,:) = squeeze(sigma(i,l,:)) + squeeze(1/Fjacb).*squeeze(F(i,k,:)).*squeeze(S(k,j,:)).*squeeze(F(l,j,:));
                        end
                    end
                end
            end
        end
        
        
%         %% Update material properties
%         function sigma = updateSigma(obj,coord,cartd0)
%             [F,Fjacb,blcg,Crcg] = obj.computeDefGradient(coord,cartd0);
%             sigma = obj.computeCauchyStress(F,Fjacb,blcg,Crcg);
%         end
        
        %% Compute deformation gradient tensor
        % Compute F, b & sigma
        function [F,Fjacb,blcg,Crcg] = computeDefGradient(obj,coord,cartd0)
            % 2D
            coord = coord(:,1:2);
            
            % Coordinate's vectorization
            x = coord(obj.connec(:,:)',:)';
            x = reshape(x,[2,obj.nnode,obj.nelem]);
            
            % Deformation gradient tensor
            F = repmat(eye(3),[1 1 obj.nelem]);
            
            for i = 1:2 % 2D
                f = zeros(2,obj.nnode,obj.nelem);
                for j = 1:2
                    for a = 1:obj.nnode
                        inc = x(i,a,:).*cartd0(j,a,:);
                        f(j,a,:) = f(j,a,:) + inc;
                    end
                    F(i,j,:) = sum(f(j,:,:));
                end
            end
            
            % Determinant vectorization (Jacobian)
            [~,Fjacb] = multinverse3x3(F);
            
            % Left-Cauchy deformation tensor
            blcg = zeros(3,3,obj.nelem);
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        blcg(i,j,:) = squeeze(blcg(i,j,:)) + (squeeze(F(i,k,:))).*(squeeze(F(j,k,:)));
                    end
                end
            end
            
            % Right-Cauchy deformation tensor
            Crcg = zeros(3,3,obj.nelem);
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        Crcg(i,j,:) = squeeze(Crcg(i,j,:)) + (squeeze(F(k,i,:))).*(squeeze(F(k,j,:)));
                    end
                end
            end
        end
          
    end
end