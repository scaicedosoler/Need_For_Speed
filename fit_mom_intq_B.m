% Need for Speed
% Caicedo-Pearce
% Moments fit
% Summer 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all;

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

%Saving runs
par.opt.save_last_run=1;

warning('off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Model version
str_mod='_intq_B';
save_str='_ipc3';

%Choose parameters to estimate
par_est='par_chib_50';

%Baseline estimation
update_baseline=1;
base_est='sol_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib_50.mat'; 

%Other calibration options
L_I_cal=0; %Change the calibrated L_I to the one measured in Compustat
pat_val_sales_win=0; %Use winzorized patent value
rb_pat_val_sales_macro=1; % patent value over sales in the data computed from macro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose period to estimate
iyear=2010; %1980; %1990; %1980; %
fyear=iyear+5;

%Load initial parameters
par=par0_fun(par);
par.tol=1e-5;

load(['data_mom' save_str '_' num2str(iyear) '_' num2str(fyear) '.mat'],'dmom')


% load('sol_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib0.mat')
% 
% par=smm.par;
% eq=smm.eq;
% 
% 
% %Modify Parameters
% 
% par.chi_b=12.1; %to match 50% pib/pi_tot
% par.chi=0.057; % To get x_l
% par.lambda=1.39; %level of growth
% par.chi_e=0.0234; % level of pat_ent
% par.alpha_e=0.41; % Relative inv_ent vs pat_ent, also changes with chi_b
% par.lambda_e=0.65; %q_ent_inc
% par.nu=0.4; %Increase concentration
% par.gamma_q=-0.0595; %q_10_90

% %Best so far
% load('min_score_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib_50.mat')

% %-----------------------------------------------------
% % EC: Aug 2024--50% additional benefits 2010
% %-----------------------------------------------------
% 
% load('sol_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib_50.mat')
% 
% par=smm.par;
% eq=smm.eq;
% 
% par.iyear=2010; 
% 
% 
% %Modify Parameters
% 
% par.chi_b=24; %to match rb_pat_val_sales
% par.chi=0.072; % To get pat_f_rb
% par.L_I=0.044; % x_l 
% par.lambda=0.72; %level of growth
% par.chi_e=0.01; % level of pat_ent
% par.alpha_e=0.3; % Relative inv_ent vs pat_ent, also changes with chi_b
% par.lambda_e=0.46; %q_ent_inc
% par.nu=0.62; %Increase concentration
% par.gamma_q=-0.085; %q_10_90

%-----------------------------------------------------
% EC: Aug 2024--50% additional benefits 2010
%-----------------------------------------------------

load('sol_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib_50.mat')


%Update baseline
if update_baseline==1

    %Load baseline
    load(base_est);

    %Save xbar of baseline as a parameter
    par.xbarb=smm.eq.xbar;

    %Save V relative to quality of baseline as a parameter
    par.V_dq_b=(smm.eq.Vbar/smm.eq.dq);

    %Save V relative to quality of baseline as a parameter
    par.inv_pat_b=(smm.par.L_I/(smm.eq.x+smm.eq.xe));

    %16. Change in Patent value over sale relative to baseline
    par.pat_val_sales_b=smm.eq.mom.pat_val_sales;

end

%Save initial year in parameters
par.iyear=iyear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load baseline
iyearb=1980; fyearb=1985;

load(['data_mom' save_str '_' num2str(iyearb) '_' num2str(fyearb) '.mat'])

dmomb=dmom;

%Load data
load(['data_mom' save_str '_' num2str(iyear) '_' num2str(fyear) '.mat'])

%Create moments relative to baseline
dmom.rb_lxi_lq=dmom.lxi_lf_cit3_res/dmomb.lxi_lf_cit3_res;
dmom.rb_xi_lq=dmom.xi_lf_cit3_res/dmomb.xi_lf_cit3_res;
dmom.rb_inv_pat=dmom.inv_pat/dmomb.inv_pat;
dmom.rb_pat_val_sales=dmom.pat_val_sales/dmomb.pat_val_sales;

dmom.pib_pi_tot=0.5; %Additional benefit is 50% of total profits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other calibration options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cj_cal=0;

if cj_cal==1
    % Other parameter calibration
    if iyear==1980
        par.L_I=0.0033; % 1980-1985 from Jones (2016)
    end
    if iyear==2010
        par.L_I=0.0066; % 2010-2015 from Jones (2016)
    end

    str_mod=[str_mod '_cj'];
end

if L_I_cal==1
    if iyear==1980
        par.L_I=0.0276; % 1980-1985 Inventors/total workers from Compustat data
    end

    if iyear==2010
        par.L_I=0.03; % 2010-2015 Inventors/total workers from Compustat data
    end

    str_mod=[str_mod '_LI'];
end

if pat_val_sales_win==1
    dmom.pat_val_sales=dmom.pat_val_sales_win;
    str_mod=[str_mod '_win'];
end

if rb_pat_val_sales_macro==1
    dmom.rb_pat_val_sales=dmom.pat_val_sales_macro/dmomb.pat_val_sales_macro;
    str_mod=[str_mod '_macro'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulated Method of Moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial parameters
par0v(1)=par.lambda;
par0v(2)=par.chi_e;
par0v(3)=par.lambda_e;
par0v(4)=par.alpha_e;
par0v(5)=par.gamma_q;
par0v(6)=par.chi;
par0v(7)=par.nu;
par0v(8)=par.L_I;

if strcmp(par_est,'par_chib_50')
    %Let chi_b change
    par0v(9)=par.chi_b;
end

if iyear==1980
    %par.alpha_q=0.085;
    par.alpha_q=0.077;
    par.alpha_x=0.31;
    % %Normalize chi_b=0
    % par.chi_b=0;
end


if iyear==1990
    par.alpha_q=0.077;
    par.alpha_x=0.3;

    %Let chi_b change
    par0v(9)=par.chi_b;
end

if iyear==2010
    par.alpha_q=0.077;
    par.alpha_x=0.31;

    %Let chi_b change
    par0v(9)=par.chi_b;
end



par.iyear=iyear;

%Parameter restrictions
eps=1e-3;
npar=length(par0v);
lb=zeros(1,npar);
ub=ones(1,npar);
lb(1:4)=[eps,eps,eps,eps];
lb(5)=-0.5;

if npar==9
    ub(9)=50;
end

ub(1:3)=5;
ub(4)=1-eps;
ub(5)=0.5;
ub(6)=10;
ub(8)=0.1;

%Save in structure
par.par_est=par_est;

%Save for string
par_str=['_' par_est];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulated Method of Moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Model str
str_mod=[str_mod '_dbase']; % '_xrb_rb_xi_lq'; %'_alphas_h_chi_b'; %'_xrb_chib0' ; %'_xrb_lxi_lq'; % '_xrb_lxi_lq'; '_xi'; '_xex' '_xrb_g' '_full'; %'_g'; 
 
%Optimization algorithm
run_alg='ps'; %'f'; %'psw'; %'sa'; % 'gs'; % 'ga'; % 

%Saving options
par.opt.save_str=[ run_alg '_' num2str(iyear) '_' num2str(fyear)  save_str str_mod par_str];

%Save minimum score
par.opt.save_score=1;

%Moments weights: 
par.nmom=17;
par.wmom=ones(1,par.nmom);

%No weights
%par.wmom(8)=0;
par.wmom(9)=0;
par.wmom(10)=0;
par.wmom(11)=0;
%par.wmom(12)=0;
par.wmom(13)=0;
par.wmom(14)=0;
par.wmom(15)=0;
par.wmom(17)=0;  

%Base year: no target relative to baseline
if iyear==1980
    par.wmom(8)=0; %No arrival rate relative to baseline
    par.wmom(11)=0; %No lxi/lq  relative to baseline
    par.wmom(13)=0; %No rb_xi_lq
    par.wmom(14)=0; %No rb_inv_pat
    par.wmom(16)=0; %No rb_pat_val_sales_macro
    par.wmom(17)=1; %50 % target
end

%Test function
scoref=@(parv)score_fun_intq_B(parv,dmom,par);
s0=scoref(par0v);

%-------------------------
% Optimization algorithms
%-------------------------

if strcmp(run_alg,'ps')
    %Pattern search: for non smooth problems
    disp('Running: Pattern Search')
    [smm.parv_sol, smm.score] =  patternsearch(scoref,par0v,[],[],[],[],lb,ub);

end

if strcmp(run_alg,'f')
    %Fmincon: for smooth problems
    disp('Running: fmincon')
    [smm.parv_sol, smm.score] =  fmincon(scoref,par0v,[],[],[],[],lb,ub);

end

if strcmp(run_alg,'sa')
    rng(8029)% For reproducibility
    %Simulated Annealing
    disp('Running: Simulated Annealing')
    [smm.parv_sol, smm.score] =  simulannealbnd(scoref,par0v,lb,ub);

end

if strcmp(run_alg,'ga')
    rng(6372)% For reproducibility
    %Genetic algorithm
    disp('Running: Genetic algorithm')
    [smm.parv_sol, smm.score] =  ga(scoref,par0v,lb,ub);

end

if strcmp(run_alg,'psw')
    rng(9284)% For reproducibility
    %Genetic algorithm
    disp('Running: Particle Swarm')
    [smm.parv_sol, smm.score] =  particleswarm(scoref,length(par0v),lb,ub);

end

if strcmp(run_alg,'gs')
    disp('Running: Global Search')
    %Global Search (finite boundaries)
    rng(3528) % For reproducibility
    gs = GlobalSearch('MaxTime',10000); %No run of local solver takes more than 3 min
    problem = createOptimProblem('fmincon','x0',par0v,...
        'objective',scoref,'lb',lb,'ub',ub);
    smm.parv_sol = run(gs,problem);
    smm.score=scoref(smm.parv_sol);

end

%Save best solution
[smm.f,smm.par,smm.eq,smm.mom]=score_mom_intq_B(smm.parv_sol,dmom,par);
smm.iyear=iyear;

%Update baseline
if update_baseline==1

    %Save xbar of baseline as a parameter
    par.xbarb=smm.eq.xbar;

    %Save V relative to quality of baseline as a parameter
    par.V_dq_b=(smm.eq.Vbar/smm.eq.dq);

    %Save V relative to quality of baseline as a parameter
    par.inv_pat_b=(smm.par.L_I/(smm.eq.x+smm.eq.xe));

    %16. Change in Patent value over sale relative to baseline
    par.pat_val_sales_b=smm.eq.mom.pat_val_sales;

end

%Saving solution
save([par.opt.dir_mat 'sol_' par.opt.save_str '.mat' ],"smm","dmom") ;




