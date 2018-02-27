classdef Physical_Problem < FEM
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = protected)
        variables
        mesh
        dim
        dof
        bc
        problemID
        nfields
        %mover a protected
        interpolation_geometry
        quadrature_geometry
%         geometry
        
        interpolation_variable
        quadrature_variable
%         geometry_variable
    end
    
    
    %% Restricted properties definition ===================================
    properties (GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
        material
        element
        physicalVars

    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID);
            obj.interpolation_geometry = Interpolation.create('mesh');
            obj.interpolation_geometry.compute(obj.mesh);
            obj.quadrature_geometry = Quadrature (obj.interpolation_geometry.geometry_type,obj.interpolation_geometry.order);
%             obj.geometry = Geometry(obj.interpolation_geometry,obj.quadrature_geometry,obj.mesh.nelem);
            
            
           
                obj.nfields=2; % nomes per stokes, s'haura de posar a algun altre lloc
                order= {'QUADRATIC','LINEAR'};
%                 order = strjoin(order);
                quadrature_variable = Quadrature (obj.mesh.geometryType,strjoin(order(1)));
                for i=1:obj.nfields
%                 quadrature_variable(i) = Quadrature (obj.mesh.geometryType,strjoin(order(i)));
                interpolation_variable(i) = Interpolation.create ('variable');
                interpolation_variable(i).compute(obj.interpolation_geometry,strjoin(order(i)));
                geometry(i) = Geometry(interpolation_variable(i),quadrature_variable,obj.mesh.nelem); 
                end
            obj.quadrature_variable = quadrature_variable;
            obj.interpolation_variable = interpolation_variable;
            obj.geometry_variable = geometry;
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);


%             obj.geometry=Geometry(obj.mesh);
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem); 
        end
        
        function preProcess(obj)
            % Create Objects

            obj.bc = BC.create(obj.mesh.ptype);
            obj.bc.compute(obj.dim.nunkn,obj.problemID,obj.interpolation_variable,obj.interpolation_geometry,obj.geometry_variable,obj.mesh.nelem);
            
            for ifield = 1:obj.nfields
            dof(ifield) = DOF(obj.interpolation_variable(ifield),obj.dim.nunkn(ifield),obj.bc,ifield);
            end
            
            obj.dof=dof;
            obj.element = Element.create(obj.mesh.ptype,obj.mesh.pdim);
            obj.physicalVars = PhysicalVariables.create(obj.mesh.ptype,obj.mesh.pdim);
            obj.solver = Solver_Dirichlet_Conditions;
        end
        
        function computeVariables(obj)
          
            obj.element.computeLHS(obj.dim,obj.mesh.nelem,obj.geometry_variable,obj.nfields,obj.material);
            obj.element.computeRHS(obj.dim,obj.mesh.nelem,obj.geometry_variable,obj.bc,obj.dof);
            
            % Assembly
            [obj.LHS,obj.RHS] = obj.Assemble(obj.element,obj.geometry_variable,obj.dim,obj.dof,obj.bc,obj.interpolation_variable,obj.nfields);
            [global_dof,global_fixnodes] = obj.assemble_dof_bc(obj.dof,obj.bc,obj.nfields);
            % Solver
            sol = obj.solver.solve(obj.LHS,obj.RHS,global_dof,global_fixnodes);
            obj.variables = obj.physicalVars.computeVars(sol,obj.dim,obj.geometry,obj.mesh.nelem,obj.dof.idx,obj.element,obj.material);
        end
        
        function postProcess(obj)
            iter = 1; % static
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,iter,results);
        end
        
        function setMatProps(obj,props)
            obj.material = obj.material.setProps(props);
        end
        function Msmooth=computeMass(obj,job)
            meshMass=obj.mesh;
            meshMass.geometryType='Triangle_Linear_Mass';
            geom=Geometry(meshMass);
            lnods=obj.mesh.connec';
            emat = zeros(geom.nnode,geom.nnode,obj.mesh.nelem);
            for igaus=1:geom.ngaus
                for inode=1:geom.nnode
                    for jnode=1:geom.nnode
                        emat(inode,jnode,:)=squeeze(emat(inode,jnode,:)) + geom.weigp(igaus)*geom.shape(inode,igaus)*geom.shape(jnode,igaus)*geom.djacob(:,igaus);
                    end
                end
            end
            
            if (job==1)
                % lumped mass matrix
                elumped = zeros(geom.nnode,obj.mesh.nelem);
                Msmooth = zeros(geom.nnode,1);
                [nproc,coeff] = nprocedure(etype,nnode);
                if (nproc==1)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
                        end
                    end
                elseif (nproc==2)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
                        end
                        elumped(inode,:)=elumped(inode,:)*coeff(inode);
                    end
                end
                for inode=1:nnode
                    Msmooth = Msmooth + sparse(lnods(inode,:),1,elumped(inode,:),npnod,1);
                end
            elseif (job==2)
                
                Msmooth = sparse(obj.mesh.npnod,obj.mesh.npnod);
                for k=1:geom.nnode
                    for l=1:geom.nnode
                        vmass = squeeze(emat(k,l,:));
                        Msmooth = Msmooth + sparse(lnods(k,:),lnods(l,:),vmass,obj.mesh.npnod,obj.mesh.npnod);
                    end
                end
                
            end
        end
        function StifMat=computeKsmooth(obj)
            StifMat=sparse(obj.mesh.npnod,obj.mesh.npnod);
            nnode=obj.geometry.nnode;
            nunkn=1;
            nstre=2;
            element_smooth=Element.create('THERMAL',obj.mesh.pdim);
            element_smooth.computeLHS(nunkn,nstre,obj.mesh.nelem,obj.geometry,obj.material);
            estiff=element_smooth.LHS;
            lnods=obj.mesh.connec';
            for a=1:nnode
                for i=1:nunkn
                    idx(nunkn*a-nunkn+i,:) = nunkn.*lnods(a,:)-nunkn+i;
                end
            end
            for k=1:nnode*nunkn
                for l=1:nnode*nunkn
                    vestiff = squeeze(estiff(k,l,:));
                    StifMat = StifMat + sparse(idx(k,:),idx(l,:),vestiff,obj.mesh.npnod,obj.mesh.npnod);
                end
            end
        end
    end
    
    %% Private methods definition =========================================
    methods (Access = protected, Static)
        function [LHS,RHS] = Assemble(element,geometry,dim,Dof,bc,interpolation,nfields)
            
            for ifield = 1:nfields
                for jfield = 1:nfields
                    
                    idx1 = Dof(ifield).idx;
                    idx2 = Dof(jfield).idx;
                    nunkn1 = dim.nunkn(ifield);
                    nnode1 = geometry(ifield).nnode;
                    nunkn2 = dim.nunkn(jfield);
                    nnode2 = geometry(jfield).nnode;
                    col = Dof(jfield).ndof;
                    row = Dof(ifield).ndof;
                    
                    % Compute LHS
                    LHS = sparse(row,col);
                    for i = 1:nnode1*nunkn1
                        for j = 1:nnode2*nunkn2
                            mat = element.LHS{ifield,jfield};
                            a = squeeze(mat(i,j,:));
                            LHS = LHS + sparse(idx1(i,:),idx2(j,:),a,row,col);
                            
                        end
                    end
                    
                    if ifield == 1 && jfield == 1
                        LHS = 1/2 * (LHS + LHS');
                    end
                    
                    LHS_global{ifield,jfield}=LHS;
                    
                    
                end
                % Compute RHS
                RHS = zeros(Dof(ifield).ndof,1);
                nnode=geometry(ifield).nnode;
                nunkn = dim.nunkn(ifield);
                idx = Dof(ifield).idx;
                for i = 1:nnode*nunkn
                    b = squeeze(element.RHS(i,1,:));
                    ind = idx(i,:);
                    RHS = RHS + sparse(ind,1,b,Dof(ifield).ndof,1);
                end
                RHS_global{ifield,1}=RHS;
            end
            LHS = cell2mat(LHS_global);
            RHS = cell2mat(RHS_global);
            
            
            
%             %Compute Global Puntual Forces
% %             if ~isempty(bc.iN)
% %                 FextPoint = zeros(dof.ndof,1);
% %                 FextPoint(bc.iN)=bc.neunodes(:,3);
% %                 RHS = RHS + FextPoint;
% %             end
%  vL = [Dof(1).vL'; Dof(2).vL' + Dof(1).ndof];
%   vR = [Dof(1).vR; Dof(2).vR + Dof(1).ndof];
%   fixnodes = [bc.fixnodes_u(:,3); bc.fixnodes_p(:,3)];
% %   fixnodes = bc.fixnodes_u(:,3);
%   x = zeros(Dof(1).ndof+Dof(2).ndof,1);
%   x(vR) = fixnodes;
%   x(vL,1) = LHS(vL,vL)\(RHS(vL) - LHS(vL,vR)*x(vR));
% 
%   [vel,p] =analytical_sol(interpolation(1).xpoints,interpolation(2).xpoints);
%   ind=1;
%   for i=1:length(interpolation(1).xpoints(:,1))
%       sol(i,:) = [interpolation(1).xpoints(i,1) interpolation(1).xpoints(i,2) x(ind) x(ind+1)];
%       ind=ind+2;
%   end
%   
%   close all
%   
% %   figure
%    pres= [x(Dof(2).vL' + Dof(1).ndof); bc.fixnodes_p(:,3)];
% %    pres= x(Dof(2).vL' + Dof(1).ndof); 
% %    plot3(interpolation(2).xpoints(Dof(2).vL,1),interpolation(2).xpoints(Dof(2).vL,2),pres','+')
%    
%    
%    
%   figure;
%   x = sol(:,1);
%   y = sol(:,2);
%   u = sol(:,3);
%   v = sol(:,4);
% %   quiver(x,y,u,v)
% 
% 
% xp = [interpolation(2).xpoints(Dof(2).vL',1); interpolation(2).xpoints(Dof(2).vR',1) ];
% yp = [interpolation(2).xpoints(Dof(2).vL',2) ;interpolation(2).xpoints(Dof(2).vR',2) ];
% % xp = interpolation(2).xpoints(:,1); 
% % yp = interpolation(2).xpoints(:,2);
% 
% tri_p = delaunay(xp,yp);
% trisurf(tri_p,xp,yp,pres)
% 
% figure
% trisurf(tri_p,xp,yp,p)
% title('P_analitica')


%     figure
%    tri= delaunay(x,y);
% %   triplot(tri,x,y,'Color',[1,1,1]*0.75)
% %   hold on; quiver(x,y,u,v,2,'k'); hold off;
%   x0=linspace(min(x),max(x),20);
%   y0=linspace(min(x),max(x),20);
% %   y0= ones(20,1)*(min(y)+max(y))/2+.;
% FlowP=TriStream(tri,x,y,u,v,x0,y0);
% PlotTriStream(FlowP,'r');
% 
% figure
% FlowP=TriStream(tri,x,y,vel(1:2:end),vel(2:2:end),x0,y0);
% PlotTriStream(FlowP,'r');
% % hold on; plot(x0,y0,'r.'); hold off;
% %   figure;
%   startx= x;
%   starty = y;
% %  streamline(x,y,u,v,startx,starty)
% 
% figure
% trisurf(tri,x,y,u)
% title('u_h')
% 
% figure
% trisurf(tri,x,y,vel(1:2:end))
% title('u_{analitica}')
% 
% figure
% trisurf(tri,x,y,v)
% title('v_h')
% 
% figure
% trisurf(tri,x,y,vel(2:2:end))
% title('v_{analitica}')
% 
% figure
% trisurf(tri,x,y,sqrt(u.^2+v.^2))
% title('modul ')
% 
% figure
% trisurf(tri,x,y,sqrt(vel(1:2:end).^2+vel(2:2:end).^2))
% title('modul_{analitica}')


        end
        
       
    end
    methods (Static)
        function [global_dof,global_fixnodes] = assemble_dof_bc(dof,bc,nfields)
            if nfields >1
                global_dof.ndof=0;
                for ifield = 1:nfields
                    global_dof.vL{ifield,1} = dof(ifield).vL' + global_dof.ndof;
                    global_dof.vR{ifield,1} = dof(ifield).vR + global_dof.ndof;
                    global_dof.ndof = global_dof.ndof + dof(ifield).ndof;
                end
                global_dof.vL = cell2mat(global_dof.vL); 
                global_dof.vR = cell2mat(global_dof.vR);
                
                global_fixnodes = [bc.fixnodes_u;bc.fixnodes_p];

            else
                global_dof.vL=dof.vL;
                global_dof.vR=dof.vR;
                global_dof.ndof = dof.ndof;
                
                global_fixnodes = bc.fixnodes;
            end
        end
    end
    
end

