classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    %% !! NEEDS REVISION !! -> should B be a class?? Or just be contained in element ??
    
    properties (GetAccess = {?Physical_Problem, ?Element_Elastic_Micro,?Element_stokes,?Element_Elastic}, SetAccess = protected)
        RHS
        LHS
    end
    
    properties (GetAccess = {?Element_Elastic,?Element_Thermal,?PhysicalVariables,?Element_Elastic_Micro,?Element_stokes}, SetAccess = {?Physical_Problem,?Element, ?Element_Elastic_Micro})
        B
    end
    
    methods (Access = ?Physical_Problem, Static)
        function element = create(ptype,pdim)
            switch ptype
                case 'ELASTIC'
                    switch pdim
                        case '2D'
                            element = Element_Elastic;
                            element.B = B2;
                        case '3D'
                            element = Element_Elastic;
                            element.B = B3;
                    end
                case 'THERMAL'
                    element = Element_Thermal;
                    element.B = B_thermal;
                    
                case 'Stokes'
                    element = Element_stokes;
                    element.B = B_stokes;
                    
                otherwise
                    error('Invalid ptype.')
            end
        end
    end
    methods (Access = {?Physical_Problem, ?Element})
        function obj = computeRHS(obj,dim,nelem,nnode,bc,dof)
%             RHSSuperficial  = obj.computeSuperficialRHS(nunkn,nelem,nnode,bc,idx);
            RHSVolumetric  = obj.computeVolumetricRHS(dim,nelem,nnode,bc,dof);
%             obj.RHS = RHSSuperficial + RHSVolumetric;
            obj.RHS = RHSVolumetric;
        end
        
        function obj = computeLHS(obj,dim,nelem,geometry_variable,nfields,material)
            for ifield=1:nfields
                for jfield=1:nfields
                    obj.LHS{ifield,jfield}= obj.computeMatrix(dim,nelem,geometry_variable(ifield),geometry_variable(jfield),material,ifield,jfield);
                end
            end
%             obj.LHS = LHS;
        end
    end
    
    methods (Abstract, Access = protected)
        r = computeSuperficialRHS(obj)
        r = computeVolumetricRHS(obj)
        r = computeMatrix(obj)
    end
end

