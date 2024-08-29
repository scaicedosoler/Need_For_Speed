% Need for Speed: Innovation Quality and the Allocation of Inventors
% Caicedo-Pearce
% Internal with scaling and additional benefit
% Main: Tables and Figures
% Summer 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear all;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-calibrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Additional benefit equal to 50% of profits

%% EC: Aug 2024--50% additional benefits
clc; clear all; close all;

load('sol_ps_1980_1985_ipc3_intq_B_macro_dbase_par_chib0.mat')

par=smm.par;
eq=smm.eq;

par.iyear=1980; 


%Modify Parameters

par.chi_b=12.1; %to match 50% pib/pi_tot
par.chi=0.057; % To get x_l
par.lambda=1.39; %level of growth
par.chi_e=0.0234; % level of pat_ent
par.alpha_e=0.41; % Relative inv_ent vs pat_ent, also changes with chi_b
par.lambda_e=0.65; %q_ent_inc
par.nu=0.4; %Increase concentration
par.gamma_q=-0.0595; %q_10_90

%par.L_I=0.035; % x_l 

%Load baseline
iyearb=1980; % 2010; %
fyearb=iyearb+5;

load(['data_mom_ipc3_' num2str(iyearb) '_' num2str(fyearb) '.mat'])

dmomb=dmom;

% Data load
load(['data_mom_ipc3_' num2str(par.iyear) '_' num2str(par.iyear+5) '.mat'])

%Create moments relative to baseline
dmom.rb_lxi_lq=dmom.lxi_lf_cit3_res/dmomb.lxi_lf_cit3_res;
dmom.rb_xi_lq=dmom.xi_lf_cit3_res/dmomb.xi_lf_cit3_res;
dmom.rb_inv_pat=dmom.inv_pat/dmomb.inv_pat;
dmom.rb_pat_val_sales=dmom.pat_val_sales_macro/dmomb.pat_val_sales_macro;


%Run equilibrium
disp('-----------------------------------------------------------')
disp('Model:Internal with scaling and additional benefit')
disp('----------------------------------------------------------')

%Equilibrium
eq1=eq_sim_fun_intq_B(par);


%Moments
mom=eq1.mom;

dmomt.g=dmom.g;
dmomt.inv_ent=dmom.inv_ent;
dmomt.pat_ent=dmom.pat_ent;
dmomt.q_ent_inc=dmom.q_ent_inc;
dmomt.inv_10=dmom.inv_10;
dmomt.pat_10=dmom.pat_10;
dmomt.q_10_90=dmom.q_10_90;
dmomt.pat_f_rb=dmom.pat_f_rb;
dmomt.x_l=dmom.x_l;
dmomt.pat_val_sales=dmom.pat_val_sales;
dmomt.pat_val_sales_win=dmom.pat_val_sales_win;
dmomt.rb_lxi_lq=dmom.rb_lxi_lq;
dmomt.rb_xi_lq=dmom.rb_xi_lq;
dmomt.rb_inv_pat=dmom.rb_inv_pat;
dmomt.rb_pat_val_sales_macro=dmom.rb_pat_val_sales;

%Display comparison
{'g',mom.g,dmom.g;
 'inv_ent',mom.inv_ent,dmom.inv_ent;
 'pat_ent',mom.pat_ent,dmom.pat_ent;
 'q_ent_inc',mom.q_ent_inc,dmom.q_ent_inc;
 'inv_10',mom.inv_10,dmom.inv_10;
 'pat_10',mom.pat_10,dmom.pat_10;
 'q_10_90',mom.q_10_90,dmom.q_10_90;
 'pat_f_rb',mom.pat_f_rb,dmom.pat_f_rb;
 'x_l',mom.x_l,dmom.x_l;
 }

%Additional benefits for profits
eq1.pib/(eq1.pi+eq1.pib+eq1.piq)



%Table of variables
%tab_fun(par,eq) % par_80 ; par_90 ; par_10

