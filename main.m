% Need for Speed: Innovation Quality and the Allocation of Inventors
% Caicedo-Pearce
% Internal with scaling and additional benefit
% Main: Tables and Figures
% Fall 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear;
warning('off')

%Figure format
par=fig_format();

%Path
addpath('./mat/')
par.opt.dir_fig='./figures/';
par.opt.dir_tab='./tables/';
par.opt.dir_mat='./mat/';

%Saving options
par.opt.save_fig=0;
par.opt.save_tab=0;

%Choose estimation version
est_str='macro_dbase_2d_alphas';%'macro_dbase_rev_alphas'; %'macro_dbase'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading parameters and moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load parameters
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_80=smm.par;dmom_80=dmom;smm_80=smm;

%load(['sol_ps_1990_1995_ipc3_intq_B_' est_str '_par_chib0.mat'])
load('sol_ps_1990_1995_ipc3_intq_B_macro_dbase_par_chib0.mat')
par_90=smm.par;dmom_90=dmom;smm_90=smm;

load(['sol_ps_2010_2015_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_10=smm.par;dmom_10=dmom;smm_10=smm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_eq=0;

if run_eq==1
    
    par=par_80; % par_80 ; par_90 ; par_10

    par.chi_b=20;

    %Run equilibrium
    disp('-----------------------------------------------------------')
    disp('Model:Internal with scaling and additional benefit')
    disp('----------------------------------------------------------')
    eq=eq_sim_fun_intq_B(par);

    %Table of variables
    tab_fun(par,eq) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table 2: Parameters & Table 4: Moments Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Load parameters
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_80=smm.par;dmom_80=dmom;smm_80=smm;eq_80=smm_80.eq;

load(['sol_ps_2010_2015_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_10=smm.par;dmom_10=dmom;smm_10=smm;eq_10=smm_10.eq;

tab_fit_cmp_intq_B(dmom_80,smm_80,dmom_10,smm_10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Validation: Untargeted moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

par_80.opt.save_fig=1;

untarget_mom_intq_B(par_80,eq_80,par_10,eq_10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table 5: Growth Decomposition; Table 6: Allocation of Inventors; Table 7: Speed and Quality Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eq0=eq_80;par0=par_80;
eq1=eq_10;par1=par_10;

g_deco_intq_B(eq_80,eq_10,par_80,par_10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 8: Counterfactual Changing Speed/Quality Labor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load parameters
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_80=smm.par;dmom_80=dmom;smm_80=smm;eq_80=smm_80.eq;

load(['sol_ps_2010_2015_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_10=smm.par;dmom_10=dmom;smm_10=smm;eq_10=smm_10.eq;

par_10.opt.save_fig=1;

cf_lab_intq_B(par_10,eq_10,par_80,eq_80)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 9: Endogenous vs. Fixed Quality (1980-1985)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par=smm.par;

% Estimate model with no endogenous quality
dmom.g=smm.eq.g; %Match estimated growth

disp('Estimation for no endogenous quality...')
smm_nq=fit_fun_mom_intq_B_nq(par,dmom);
smm.par.lambda_nq=smm_nq.par.lambda;

smm.par.save_str=['_' smm.par.par_est '_1980_1985'];
smm.par.opt.save_fig=0;

%Counterfactuals
disp('---------------')
disp('Counterfactual')
disp('---------------')
cf_nq_fun_intq_B(smm)

% %---------------------------
% %For 2010-2015 calibration
% %---------------------------
% load(['sol_ps_2010_2015_ipc3_intq_B_macro_dbase_par_chib0.mat'])
% par=smm.par;
% 
% % Estimate model with no endogenous quality
% dmom.g=smm.eq.g; %Match estimated growth
% 
% disp('Estimation for no endogenous quality...')
% smm_nq=fit_fun_mom_intq_B_nq(par,dmom);
% smm.par.lambda_nq=smm_nq.par.lambda;
% 
% smm.par.save_str=['_base_' smm.par.par_est '_2010_2015'];
% smm.par.opt.save_fig=0;
% 
% %Counterfactuals
% disp('---------------------------------------------')
% disp('Counterfactual: No Endogneous Quality')
% disp('---------------------------------------------')
% cf_nq_fun_intq_B(smm)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 10: Robustness Initial Additional Benefit 50% of Profits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc,clear all, close all;

est_str='macro_dbase' ;%'macro_dbase_rev_alphas'; %

% Load parameters
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib_50.mat'])
par_50_80=smm.par;dmom_50_80=dmom;smm_50_80=smm;eq_50_80=smm.eq;

load(['sol_ps_2010_2015_ipc3_intq_B_' est_str '_par_chib_50.mat'])
par_50_10=smm.par;dmom_50_10=dmom;smm_50_10=smm;eq_50_10=smm.eq;

tab_fit_cmp_intq_B(dmom_50_80,smm_50_80,dmom_50_10,smm_50_10)


par_50_10.opt.save_fig=1;
par_50_10.save_str='_50';
cf_lab_intq_B(par_50_10,eq_50_10,par_50_80,eq_50_80)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 11: Counterfactual changing alpha_x and alpha_q 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc;clear ; close all

% Load parameters
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_80=smm.par;dmom_80=dmom;smm_80=smm;eq_80=smm_80.eq;

load(['sol_ps_2010_2015_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_10=smm.par;dmom_10=dmom;smm_10=smm;eq_10=smm_10.eq;

par_10.opt.save_fig=1;
cf_alpha_q_alpha_x(par_10,eq_10,par_80,eq_80)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 13: Counterfactual External Effects of Speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc; clear all; close all

% load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
% %load(['sol_ps_2010_2015_ipc3_intq_B_macro_dbase_par_chib0.mat'])
% par=smm.par;
% 
% %Save option
% par.opt.save_fig=0;

% %Computes the counterfactual of parameters depending on x
% disp('---------------------------------------------')
% disp('Counterfactual: Parameters depending on x')
% disp('---------------------------------------------')
% cf_par_x_fun_intq_B(par)

% %Computes the counterfactual of parameters depending on x
% disp('---------------------------------------------')
% disp('Counterfactual: Parameters depending on x, version 2')
% disp('---------------------------------------------')
% cf_par_x_fun_v2_intq_B(par)


%Computes the counterfactual of parameters depending on x
disp('---------------------------------------------')
disp('Counterfactual: Quality depending on x, growth effects')
disp('---------------------------------------------')
%Computes counterfactual growth for parameters
% Load parameters
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_80=smm.par;dmom_80=dmom;smm_80=smm;eq_80=smm_80.eq;

load(['sol_ps_2010_2015_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_10=smm.par;dmom_10=dmom;smm_10=smm;eq_10=smm_10.eq;

par_10.opt.save_fig=1;
cf_lab_x_intq_B(par_10,eq_10,par_80,eq_80)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table C.1: Parameters & Table C.2: Moments Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['sol_ps_1990_1995_ipc3_intq_B_' est_str '_par_chib0.mat'])
par_90=smm.par;dmom_90=dmom;smm_90=smm;

tab_fit_cmp_intq_B(dmom_80,smm_80,dmom_90,smm_90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure C.2: External Effects of Speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
%load(['sol_ps_2010_2015_ipc3_intq_B_macro_dbase_par_chib0.mat'])
par=smm.par;

%Save option
par.opt.save_fig=0;

% %Computes the counterfactual of parameters depending on x
% disp('---------------------------------------------')
% disp('Counterfactual: Parameters depending on x')
% disp('---------------------------------------------')
% cf_par_x_fun_intq_B(par)

%Computes the counterfactual of parameters depending on x
disp('---------------------------------------------')
disp('Counterfactual: Parameters depending on x, version 2')
disp('---------------------------------------------')
cf_par_x_fun_v2_intq_B(par)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Work Zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional Exercises: Counterfactual: Optimal Growth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For 1980-1985 calibration
load(['sol_ps_1980_1985_ipc3_intq_B_' est_str '_par_chib0.mat'])
par=smm.par;

%Plots relative to optimal growth
par.opt.save_fig=0;
par.opt.save_str2=[est_str '_par_chib0_'];
cf_opt_g_intq_B(par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_reg==1

    clc; clear all; close all;

    % Load parameters
    load(['sol_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib0.mat'])
    par_80=smm.par;dmom_80=dmom;smm_80=smm;eq_80=smm_80.eq;

    load(['sol_ps_2010_2015_ipc3_intq_B_macro_dbase_par_chib0.mat'])
    par_10=smm.par;dmom_10=dmom;smm_10=smm;eq_10=smm_10.eq;

    %Run equilibrium
    disp('-----------------------------------------------------------')
    disp('Model:Internal with scaling and additional benefit')
    disp('----------------------------------------------------------')

    Lg=0.01:0.005:0.1;
    ng=length(Lg);

    par=par_80;
    for ig=1:ng
        par.L_I=Lg(ig);
        eq=eq_sim_fun_intq_B(par);
        x(ig)=eq.x;
        lx(ig)=eq.lx;
        dq(ig)=eq.dq;
        lq(ig)=eq.lq;
    end

    %Coefficients
    dlogx_dlogL=diff(log(x))/diff(log(Lg));
    dlogdq_dlogL=diff(log(dq))/diff(log(Lg));


    mdl_x = fitlm(diff(log(Lg)),diff(log(x)));
    mdl_q = fitlm(diff(log(Lg)),diff(log(dq)));


% Unsure why regressions do not work! Check: S Oct 20, 2024
    % X_L=[ones(ng,1), log(Lg)'];
    % X_lx=[ones(ng,1), log(lx)'];
    % X_lq=[ones(ng,1), log(lq)'];
    % 
    % log_x=log(x)';
    % log_dq=log(dq)';
    % 
    % %Regressions
    % bx_L=log_x\X_L;
    % bx_lx=log_x\X_lx;
    % bdq_L=log_dq\X_L;
    % bdq_lq=log_dq\X_lq;

end


