classdef TopOpt_Problem_Chomog < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Chomog(settings)
            obj@TopOpt_Problem(settings);
            switch settings.ptype
                case 'Chomog_alphabeta_st_Volume'
                    obj.cost=ShFunc_Chomog_alphabeta(settings);
                    obj.constraint=ShFunc_Volume(settings);
                case 'Chomog_fraction_st_Volume'
                    obj.cost=ShFunc_Chomog_fraction(settings);
                    obj.constraint=ShFunc_Volume(settings);
                case 'ChomogLamPerimeter_alphabeta_st_Volume'
                    obj.cost=ShFunc_ChomogLamPerimeter_alphabeta(settings);
                    obj.constraint=ShFunc_Volume(settings);
                case 'ChomogLamPerimeter_fraction_st_Volume'
                    obj.cost=ShFunc_ChomogLamPerimeter_fraction(settings);
                    obj.constraint=ShFunc_Volume(settings);
                case 'ChomogLamPerimeter_st_VolumeEnforceCh_inf'
                    obj.cost=ShFunc_ChomogLamPerimeter_CC(settings);
                    obj.constraint=ShFunc_VolumeEnforceCh_inf(settings);
                otherwise
                    fprintf('Micro problem not implemented');
            end
            
        end
    end
end
