classdef ShFunc_VolumeEnforceCh_inf< Shape_Functional
    properties
        volume
        enforceCh
    end
    methods
        function obj=ShFunc_VolumeEnforceCh_inf(settings)
            %        obj@Shape_Functional(settings);
            obj.volume=ShFunc_Volume(settings);
            for i=1:6
                obj.enforceCh{i}=ShFunc_Chomog_EnforceCh_CCstar(settings,i);
                if isequal(i,5) || isequal(i,4)
                    obj.enforceCh{i}.setEpsilon(0);
                end
            end
        end
        
        function computef(obj, x, physicalProblem, interpolation,filter)
            obj.volume.target_parameters=obj.target_parameters;
            obj.volume.computef(x, physicalProblem, interpolation,filter);
            obj.value(1)=obj.volume.value;
            obj.gradient(:,1)=obj.volume.gradient;            
            for i=1:6
                obj.enforceCh{i}.target_parameters=obj.target_parameters;
                obj.enforceCh{i}.computef(x,physicalProblem,interpolation,filter);
                obj.value(i+1,1)=obj.enforceCh{i}.value;
                obj.gradient(:,i+1)=obj.enforceCh{i}.gradient;
            end
        end
    end
end