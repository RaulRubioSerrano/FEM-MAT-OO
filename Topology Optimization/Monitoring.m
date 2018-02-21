classdef Monitoring < handle
    %Monitoring Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ON
        figures
        monitor
        nfigs
        ncost
        nconstraint
        nstop
        stop_names
    end
    
    methods
        function obj = Monitoring(settings,ON)
            obj.ON = ON;
            if obj.ON
                obj.getStopVarsNames(settings.optimizer);
                
                obj.figures = {};
                obj.createFigure('Cost');
                
                obj.ncost = length(settings.cost);
                for i = 1:obj.ncost
                    if isempty(settings.weights)
                        obj.createFigure([obj.setCase(settings.cost{i}) ' (wt. 1.0)']);
                    else
                        obj.createFigure([obj.setCase(settings.cost{i}) sprintf(' (wt. %.2f)',settings.weights(i))]);
                    end
                end
                
                obj.nconstraint = length(settings.constraint);
                for i = 1:obj.nconstraint
                    obj.createFigure(['Cstr ' num2str(i) ': ' obj.setCase(settings.constraint{i})]);
                    obj.createFigure(['\lambda' obj.subindex(obj.setCase(settings.constraint{i}))]);
                end
                
                for i = 1:obj.nstop
                    obj.createFigure(['Conv Criteria ' num2str(i) ': ' obj.stop_names{i}]);
                end
                
                obj.nfigs = length(obj.figures);
                if obj.nfigs <= 4
                    nrows = 1;
                else
                    nrows = 2;
                end
                ncols = round(obj.nfigs/nrows);
                
                % Create obj.figures
                obj.monitor = figure;
                
                % WARNING!! JavaFrame will be obsolete in future Matlab releases
                warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
                drawnow; set(get(obj.monitor,'JavaFrame'),'Maximized',1);
                
                for i = 1:obj.nfigs
                    obj.figures{i}.style = obj.subplot_tight(nrows,ncols,i,[0.06 0.03]);
                    switch obj.figures{i}.chart_type
                        case 'plot'
                            obj.figures{i}.handle = plot(0,0);
                            title(obj.figures{i}.title);
                            grid on
                        case 'bar'
                            obj.figures{i}.handle = bar(0,0);
                            title(obj.figures{i}.title);
                            grid on
                    end
                end
            end
        end
        
        function display(obj,iteration,cost,constraint,stop_vars)
            if obj.ON
                set(obj.monitor,'NumberTitle','off','Name',sprintf('Monitoring - Iteration: %.0f',iteration))
                obj.figures{1} = obj.updateFigure(obj.figures{1},iteration,cost.value);
                for i = 1:obj.ncost
                    k = i+1;
                    obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,cost.ShapeFuncs{i}.value);
                end
                for i = 1:obj.nconstraint
                    k = i*2+obj.ncost;
                    obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,constraint.ShapeFuncs{i}.value);
                    k = k+1;
                    obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,constraint.lambda(i));
                end
                for i = 1:obj.nstop
                    k = i+1+obj.ncost+2*obj.nconstraint;
                    if contains(obj.figures{k}.title,'outit')
                        obj.updateFigureIT(obj.figures{k},iteration,stop_vars(i,:));
                    else
                        obj.figures{k} = obj.updateFigure(obj.figures{k},iteration,stop_vars(i,1));
                    end
                end
            else
                return
            end
        end
    end
    methods (Access = private)
        function obj = createFigure(obj,TITLE)
            obj.figures{end+1}.title = TITLE;
            obj.figures{end}.iteration = [];
            obj.figures{end}.variable = [];
            if contains(TITLE,'kappa') || contains(TITLE,'outit')
                obj.figures{end}.chart_type = 'bar';
            else
                obj.figures{end}.chart_type = 'plot';
            end
        end
        
        function obj = getStopVarsNames(obj,optimizer)
            switch optimizer
                case 'SLERP'
                    obj.stop_names = {'\Deltacost';'Norm L2';'\kappa'};
                case 'PROJECTED GRADIENT'
                    obj.stop_names = {'\Deltacost';'Norm L2';'\kappa'};
                case 'MMA'
                    obj.stop_names = {'kktnorm';'outit'};
                case 'IPOPT'
                    obj.stop_names = {};
            end
            obj.nstop = length(obj.stop_names);
        end
    end
    methods (Access = private, Static)
        function figures = updateFigure(figures,iteration,variable)
            figures.iteration(end+1) = iteration;
            figures.variable(end+1) = variable;
            set(figures.handle,'XData',figures.iteration,'YData',figures.variable);
            set(figures.style,'XLim',[0 figures.iteration(end)])
            drawnow
        end
        
        function updateFigureIT(figures,iteration,variable)
            set(figures.handle,'XData',iteration,'YData',variable(1));
            set(figures.style,'XLim',[iteration-1 iteration+1],'YLim',[0 variable(2)])
            drawnow
        end
        
        function str = setCase(str)
            str = [upper(str(1)) lower(str(2:end))];
        end
        
        function str = subindex(str_)
            str = [];
            for i = 1:length(str_)
                str(end+1:end+2) = ['_' str_(i)];
            end
        end
        
        function vargout = subplot_tight(m, n, p, margins, varargin)
            %% subplot_tight
            % A subplot function substitude with margins user tunabble parameter.
            %
            %% Syntax
            %  h=subplot_tight(m, n, p);
            %  h=subplot_tight(m, n, p, margins);
            %  h=subplot_tight(m, n, p, margins, subplotArgs...);
            %
            %% Description
            % Our goal is to grant the user the ability to define the margins between neighbouring
            %  subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the
            %  margins between subplots can reach 40% of figure area, which is pretty lavish. While at
            %  the begining the function was implememnted as wrapper function for Matlab function
            %  subplot, it was modified due to axes del;etion resulting from what Matlab subplot
            %  detected as overlapping. Therefore, the current implmenetation makes no use of Matlab
            %  subplot function, using axes instead. This can be problematic, as axis and subplot
            %  parameters are quie different. Set isWrapper to "True" to return to wrapper mode, which
            %  fully supports subplot format.
            %
            %% Input arguments (defaults exist):
            %   margins- two elements vector [vertical,horizontal] defining the margins between
            %        neighbouring axes. Default value is 0.04
            %
            %% Output arguments
            %   same as subplot- none, or axes handle according to function call.
            %
            %% Issues & Comments
            %  - Note that if additional elements are used in order to be passed to subplot, margins
            %     parameter must be defined. For default margins value use empty element- [].
            %  -
            %
            %% Example
            % close all;
            % img=imread('peppers.png');
            % figSubplotH=figure('Name', 'subplot');
            % figSubplotTightH=figure('Name', 'subplot_tight');
            % nElems=17;
            % subplotRows=ceil(sqrt(nElems)-1);
            % subplotRows=max(1, subplotRows);
            % subplotCols=ceil(nElems/subplotRows);
            % for iElem=1:nElems
            %    figure(figSubplotH);
            %    subplot(subplotRows, subplotCols, iElem);
            %    imshow(img);
            %    figure(figSubplotTightH);
            %    subplot_tight(subplotRows, subplotCols, iElem, [0.0001]);
            %    imshow(img);
            % end
            %
            %% See also
            %  - subplot
            %
            %% Revision history
            % First version: Nikolay S. 2011-03-29.
            % Last update:   Nikolay S. 2012-05-24.
            %
            % *List of Changes:*
            % 2012-05-24
            %  Non wrapping mode (based on axes command) added, to deal with an issue of disappearing
            %     subplots occuring with massive axes.
            
            %% Default params
            isWrapper=false;
            if (nargin<4) || isempty(margins)
                margins=[0.04,0.04]; % default margins value- 4% of figure
            end
            if length(margins)==1
                margins(2)=margins;
            end
            
            %note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
            [subplot_col,subplot_row]=ind2sub([n,m],p);
            
            
            height=(1-(m+1)*margins(1))/m; % single subplot height
            width=(1-(n+1)*margins(2))/n;  % single subplot width
            
            % note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
            subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot
            subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot
            
            merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
            merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width
            
            merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
            merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
            pos=[merged_left, merged_bottom, merged_width, merged_height];
            
            
            if isWrapper
                h=subplot(m, n, p, varargin{:}, 'Units', 'Normalized', 'Position', pos);
            else
                h=axes('Position', pos, varargin{:});
            end
            
            if nargout==1
                vargout=h;
            end
        end
    end
    
end

