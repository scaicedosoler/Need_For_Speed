%Creates initial parameter structure
function par=par0_fun(par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fixed parameters
par.Lp=1; par.rho=0.02; par.beta=0.1090; par.varrho=0.9; par.qbar=1; 

%Maximum number of product lines
par.Nmax=15;

%Benchmark--Int--q-- Estimated for Int
par.L_I=0.01;
par.alpha_x=0.3; par.chi=1;  
par.lambda=0.124; par.alpha_q=0.1;
par.lambda_e=0.057; par.alpha_e=0.018; par.gamma_e=1;

% Creative destruction model: tau=Ce*(chi_e+chi_e0*Exp(-chi_e1*x))
par.chi_e=0.027; 
par.chi_e0=0.01; %0;%
par.chi_e1=1;

%Speed and quality scaling parameters
c= 1;%1.1; %1;% -9; %0%
par.gamma_q=-c*par.alpha_q; %-par.alpha_q; %0; %
par.gamma_x=1-par.alpha_x-par.alpha_q-par.gamma_q; %To make xQ linear in q (when chi_e0=0)

% Algorithm parameters and options
%Options for solvers
par.alg.opt_fz = optimset('Display','off');
par.alg.opt_fs = optimoptions('fsolve','Display','off');
par.alg.opt_fmincon0=optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
par.alg.opt_fmincon=optimoptions(@fmincon,'Algorithm','interior-point','Display','off','OutputFcn',@stopfun);

par.alg.maxIter=10000;
par.alg.crit1=1e-5;

%Equilibrium and solver tolerance
par.alg.tol=1e-4;
par.alg.Vtol=1e-4;
par.alg.Vmax_it=5000; %Maximum iterations for value function
par.alg.Vopts = optimset('TolX', 1e-5); %Optimization for bisection method

%Solver options
par.opt.fmincon=1;

%Save starting point for computing the equilibrium
par.opt.save_var0=0;

% Simulation parameters
par.sim.years=20000; %Total years
par.sim.dt=1/2; %1/2; % Discretization
par.sim.Tmax=ceil(par.sim.years/par.sim.dt); %Simulation points
par.sim.nprod=10000; %Number of product lines (to start)
par.sim.seed=22; 
par.sim.tol=1e-4;

end