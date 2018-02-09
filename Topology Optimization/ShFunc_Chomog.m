classdef ShFunc_Chomog < Shape_Functional
    properties (Access = protected)
        h_C_0
        Chomog
        tstress
        tstrain
        Chomog_Derivatives
    end
    methods
        function obj=ShFunc_Chomog(settings)
           obj.target_parameters = settings.target_parameters;            
        end
    end
    methods (Access = protected)
        function compute_Chomog_Derivatives(obj,nstre,nelem,ngaus,x,interpolation,filter)
            rho=filter.getP0fromP1(x);
            matProps=interpolation.computeMatProp(rho);
            %           mass=filter.Msmooth;
            %             obj.tstrain = permute(obj.tstrain,[2 3 4]);
            %             obj.tstress = permute(obj.tstrain,);
            obj.Chomog_Derivatives = zeros(nstre,nstre,ngaus,nelem);
            for istreChomog = 1:nstre
                for jstreChomog = 1:nstre
                    for igaus=1:ngaus
                        for istre=1:nstre
                            for jstre = 1:nstre
                                obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:) = ...
                                    squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,igaus,:)) + ...
                                    (squeeze(obj.tstrain(istreChomog,igaus,istre,:))...
                                    .*squeeze(matProps.dC(istre,jstre,:,igaus))...
                                    .*squeeze(obj.tstrain(jstreChomog,igaus,jstre,:)));
                            end
                        end
                    end
                    %                     C_D = filter.getP1fromP0(squeeze(obj.Chomog_Derivatives(istreChomog,jstreChomog,:,:)));
                    %                     obj.Chomog_Derivatives(istreChomog,jstreChomog,:,:) = mass*C_D;
                end
            end
            
        end
        
        function r = projection_Chomog(obj,inv_matCh,alpha,beta)
            weights = alpha*beta';
            r = sum(sum(weights.*inv_matCh));
        end
        
        function r = derivative_projection_Chomog(obj,inv_matCh,alpha,beta,Chomog_Derivatives,nelem,ngaus,nstre)
            weights = alpha*beta';
            weights_inv = inv_matCh*weights*inv_matCh;
            DtC1 = zeros(ngaus,nelem);
            DtC = zeros(ngaus,nelem);
            for igaus=1:ngaus
                for i=1:nstre
                    for j=1:nstre
                        DtC1(igaus,:) = squeeze(Chomog_Derivatives(i,j,igaus,:));
                        DtC(igaus,:) = DtC(igaus,:)- weights_inv(i,j)*DtC1(igaus,:);
                    end
                end
            end
            r = DtC;
        end
        
        function setPhysicalData(obj,variables)
            obj.Chomog = variables.Chomog;
            obj.tstrain = variables.tstrain;
            obj.tstress = variables.tstress;            
        end
        
         function r = compute_Ch_star(obj,TOL)
            E_plus = TOL.E_plus;
            E_minus = TOL.E_minus;
            nu_plus = TOL.nu_plus;
            nu_minus = TOL.nu_minus;
            
            C_Cstar_case = 'Seba';
            
            switch C_Cstar_case
                case 'negative_poisson'
                    kappa_f = @(E,nu) E/2*(1-nu);
                    mu_f = @(E,nu) E/2*(1-nu);
                    
                    k_plus = kappa_f(E_plus,nu_plus);
                    mu_plus = mu_f(E_plus,nu_plus);
                    
                    k_minus = kappa_f(E_minus,nu_minus);
                    mu_minus = mu_f(E_minus,nu_minus);
                    
                    
                    nu = @(k,mu) (k-mu)/(k+mu);
                    E = @(k,mu) (4*k*mu)/(k+mu);
                    C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
                    
                    kappa_nu_min = k_minus;
                    mu_nu_min = mu_plus;
                    
                    nu_min = nu(kappa_nu_min,mu_nu_min);
                    E_nu_min = E(kappa_nu_min,mu_nu_min);
                    C_nu_min = C(E_nu_min,nu_min);
                    r = C_nu_min;
                    
                case 'nu_0_6' %From Sigmund Thesis% rho = 0.38
                    C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
                    nu = -0.6;
                    E = (1-nu*nu)*0.04;
                    r = C(E,nu);
                    
                case 'Seba' % Es=0.08; nus=-0.25
                    
                    r = [0.0853    -0.0213       0;
                        -0.0213    0.0853       0;
                        0         0    0.0533];
                    
                case 'nu_0_8' %From Sigmund Thesis% rho = 0.25
                    C = @(E,nu) (E/(1-nu*nu)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]);
                    nu = -0.8;
                    E = (1-nu*nu)*0.02;
                    r = C(E,nu);
            end
        end
    end
end







