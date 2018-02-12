clc
clear variables;close all;

addpath(genpath('./FEM'));
addpath(genpath('./Topology Optimization'));
% addpath(genpath(fullfile('.','FEM')));
% addpath(genpath(fullfile('.','Topology Optimization')));
%% test
%run('test.m');
clear variables;
%% settings
% settings.ptype='MACRO';
% settings.filename='TOPOPT_TEST';
% settings.filename='topopt_quad';
% settings.filename='GrippingNew';

settings.ptype='MICRO';
settings.filename='RVE_Square_Triangle';
% settings.filename='RVE_Square_Triangle_Fine';


settings.plotting=true;
settings.printing=true;
settings.maxiter = 2000;


settings.method='SIMPALL';
%settings.method='SIMP_P3';
% settings.method='SIMP_Adaptative';

settings.material='ISOTROPIC';
%settings.initial_case='full';
settings.initial_case='circle';
% settings.initial_case='horizontal';
 %settings.initial_case='square';
% settings.initial_case='feasible';
% settings.initial_case='rand';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings.ptype='Chomog_alphabeta_st_Volume';
% settings.ptype='Chomog_fraction_st_Volume';
% settings.ptype='Chomog_CC_st_Volume';
% settings.ptype='ChomogLamPerimeter_alphabeta_st_Volume';
% settings.ptype='ChomogLamPerimeter_fraction_st_Volume';
% settings.ptype='ChomogLamPerimeter_CC_st_Volume';
% settings.ptype='ChomogLamPerimeter_CC_st_VolumeEnforceCh_inf';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings.cost={'chomog_CC';'perimeter'};%'chomog_xxxxx';'compliance';'perimeter';
% settings.multipliers=[]; %all 1
settings.multipliers=[1 0.1]; %compl+lambda*perimeter
settings.constraint={'volume';'enforceCh_CCstar'};

settings.optimizer='SLERP';
%settings.optimizer='PROJECTED GRADIENT';settings.kappaMultiplier=1;
% settings.optimizer='MMA';
% settings.optimizer='IPOPT';

settings.filter='P1';%'PDE';

settings.TOL.rho_plus=1;
settings.TOL.rho_minus=0;
settings.TOL.E_plus=1;
settings.TOL.E_minus=1e-3;
settings.TOL.nu_plus=1/3;
settings.TOL.nu_minus=1/3;


settings.target_parameters.Vfrac=0.5;
settings.target_parameters.optimality_tol=1e-3;
settings.target_parameters.constr_tol=1e-3;
settings.target_parameters.Perimeter_target=5;
settings.perimeter.optimizer=settings.optimizer;
settings.target_parameters.epsilon_isotropy = 1e-3;%%%%%%

% if strcmp(problembsc.enforceCh_type,'isotropy')
%     target_parameters.epsilon_isotropy_ini = 1e-3;
%     target_parameters.epsilon_isotropy_final = 1e-5;
% end


settings.nsteps=1;
settings.Vfrac_final=settings.target_parameters.Vfrac;
settings.optimality_final=settings.target_parameters.optimality_tol;
settings.constr_final=settings.target_parameters.constr_tol;
settings.epsilon_isotropy_final=settings.target_parameters.epsilon_isotropy;
settings.Vfrac_initial=1;

settings.optimality_initial=1e-3;
settings.constr_initial=1e-3;
settings.epsilon_isotropy_initial=1e-1;

settings.micro.alpha =[1 0 0]';%[1 0 0]'
settings.micro.beta =[0 -1 0]';%[0 -1 0]'
settings.micro.selectiveC_Cstar = [1,1,1;
                            1,1,1;
                             1,1,1];

%% main

tic
test=TopOpt_Problem(settings);
% test.checkDerivative;
% toc
test.preProcess;
test.computeVariables;
toc

test.postProcess;


