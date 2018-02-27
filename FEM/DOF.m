classdef DOF
    %DOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem, ?Element, ?Solver}, SetAccess = private)
        
    end
    properties (GetAccess = {?Physical_Problem, ?Element}, SetAccess = private)
        idx
    end
    properties (GetAccess = public, SetAccess = private)
        vR
        ndof
        vL
    end
    
    methods (Access = {?Physical_Problem,?Filter})
        % Constructor
        function obj = DOF(interpolation,nunkn,bc,ifield)
            % Compute idx
            fixnodes = obj.select_fixnodes(bc,ifield);
            npnod = length(interpolation.xpoints(:,1));
            nnode = interpolation.isoparametric.nnode;
            connec = interpolation.T;
%             nunkn = nunkn.u;
            lnods = connec';
            for i = 1:nnode
                for j = 1:nunkn
                    obj.idx(nunkn*(i-1)+j,:) = nunkn.*lnods(i,:) - nunkn + j;
                end
            end
            obj.ndof = nunkn*npnod;
            
            % *************************************************************
            if (size(fixnodes,1)>0)
                vR = (fixnodes(:,1)-1)*nunkn + fixnodes(:,2);  % Finds the equation number
                vL = setdiff (1:obj.ndof, vR);
            else
                vL = (1:obj.ndof);
                vR = [];
            end
            
            % *************************************************************
%             switch linearTriangle.type
%                 case {'TRIANGLE','QUADRILATERAL'}
%                     if (size(fixnodes,1)>0)
%                         vR = (fixnodes(:,1)-1)*nunkn + fixnodes(:,2);  % Finds the equation number
%                         vL = setdiff (1:obj.ndof, vR);
%                     else
%                         vL = (1:obj.ndof);
%                         vR = [];
%                     end
%                 case 'LINEAR_TRIANGLE_MIX'
%                     if (size(fixnodes,1)>0)
%                         vR = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  % Finds the equation number
%                         vL = setdiff (1:obj.ndof, vR);
%                     else
%                         vL = (1:obj.ndof);
%                         vR = [];
%                     end
%                 case 'LINEAR_TRIANGLE_MIX_COUPLED'
%                     vR = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  % Finds the equation number
%                     vL = setdiff (1:obj.ndof, vR);
%                     
%                 case 'HEXAHEDRON'
%                     if (size(fixnodes,1)>0)
%                         vR = (fixnodes(:,1)-1)*3 + fixnodes(:,2);  % Finds the equation number
%                         vL = setdiff (1:obj.ndof, vR);
%                     else
%                         vL = (1:obj.ndof);
%                         vR = [];
%                     end
%                 otherwise
%                     error('No existe es tipo de elemento o no ha sido implementado')
%             end
            obj.vR = vR;
            obj.vL = vL;
        end
        
    end
    methods (Static)
        function fixnodes = select_fixnodes(bc,ifield)
            if ifield == 1
                fixnodes = bc.fixnodes_u;
            else
                fixnodes = bc.fixnodes_p;
            end
        end
    end
    
end

