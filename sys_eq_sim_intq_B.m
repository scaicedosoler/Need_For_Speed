%System of equations with internal innovation and additional benefit

function [f,eq]=sys_eq_sim_intq_B(var,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Unpack variables var=(lx,lq)
eq.lx=var(1);
eq.lq=var(2);

%Mass of product lines
eq.Phi=1; 

% Minimum quality
eq.qmin=0; 

%Benefit
eq.B=par.chi_b;

%Finals good production
eq.betatilde=par.beta^par.beta*(1-par.beta)^(1-2*par.beta);
eq.lp=par.beta*par.Lp/(par.beta+eq.Phi*(1-par.beta)^2);
eq.pi=eq.lp*(1-par.beta)*eq.betatilde;

eq.Y=(eq.betatilde/par.beta)*par.qbar*eq.lp;
eq.wp=eq.betatilde*par.qbar;
eq.pj=eq.wp/((1-par.beta)*par.qbar);
eq.kj=eq.lp*((1-par.beta)*par.qbar/eq.wp).^(1/par.beta);
eq.yj=eq.kj*eq.pj;

%Innovations
if par.L_I-eq.lx-eq.lq>1e-8 %Condition for le>0

    %1.Arrival rate
    eq.x=par.chi*eq.lx^par.alpha_x;

    %2.Quality
    eq.Q=par.lambda*eq.lq^par.alpha_q;

    %3.Entrants
    eq.le=par.L_I-eq.lx-eq.lq;

    %4 & 5.Creative destruction
    eq.xe=par.chi_e*eq.le^par.alpha_e;
    eq.tau=eq.xe;

    %6.Quality
    eq.dq=par.lambda*eq.lq^par.alpha_q;
    eq.dqe=par.lambda_e;

    %7.Growth
    eq.g=eq.xe*eq.dqe+eq.x*eq.dq;

    %7.Interest rate
    eq.r=eq.g+par.rho;

    %8.A 
    eq.Clq=eq.lq;
    eq.Clx=eq.lx;
    eq.pi_b=par.chi*par.chi_b*eq.Clx^(par.alpha_x)*(1-par.alpha_x);

    eq.A=(eq.pi+eq.pi_b)/((eq.r+eq.tau)-par.chi*par.lambda*eq.Clx^(par.alpha_x)*...
        eq.Clq^(par.alpha_q)*(1-par.alpha_x-par.alpha_q));

    %9. Wages
    eq.qw=eq.le^(1-par.alpha_e)/(par.alpha_e*par.chi_e*eq.A*(1+par.lambda_e));
    eq.wq=1/eq.qw;

    %--------------------------
    % Equations: (lq,lx)
    %--------------------------
    f(1)=eq.lx^(1-par.alpha_x)-par.alpha_x*eq.qw*(par.chi*par.chi_b+ ...
        (par.chi*eq.A*par.lambda)^(1/(1-par.alpha_q))*(par.alpha_q*eq.qw)^(par.alpha_q/(1-par.alpha_q))...
        *eq.lx^(par.alpha_x*par.alpha_q/(1-par.alpha_q)));
    f(2)=eq.lq-(par.alpha_q* par.chi*eq.A* par.lambda*eq.qw)^(1/(1-par.alpha_q))*...
        eq.lx^(par.alpha_x/(1-par.alpha_q));

end

if par.L_I-eq.lx-eq.lq<1e-8 %Condition for le~0

    %1.Arrival rate
    eq.x=par.chi*eq.lx^par.alpha_x;

    %2.Quality
    eq.Q=par.lambda*eq.lq^par.alpha_q;

    %3.Entrants
    eq.le=0;

    %4 & 5.Creative destruction
    eq.xe=par.chi_e*eq.le^par.alpha_e;
    eq.tau=eq.xe;

    %6.Quality
    eq.dq=par.lambda*eq.lq^par.alpha_q;
    eq.dqe=par.lambda_e;

    %7.Growth
    eq.g=eq.xe*eq.dqe+eq.x*eq.dq;

    %7.Interest rate
    eq.r=eq.g+par.rho;

    %8.A 
    eq.Clq=eq.lq;
    eq.Clx=eq.lx;
    eq.pi_b=par.chi*par.chi_b*eq.Clx^(par.alpha_x)*(1-par.alpha_x);

    eq.A=(eq.pi+eq.pi_b)/((eq.r+eq.tau)-par.chi*par.lambda*eq.Clx^(par.alpha_x)*...
        eq.Clq^(par.alpha_q)*(1-par.alpha_x-par.alpha_q));

    %9. Wages
    eq.qw=eq.lq^(1-par.alpha_q)/(par.alpha_q* par.chi*eq.A* par.lambda*eq.lx^(par.alpha_x));
    eq.wq=1/eq.qw;

    %--------------------------
    % Equations: (lq,lx)
    %--------------------------

    f(1)=eq.lx^(1-par.alpha_x)-par.alpha_x*eq.qw*(par.chi*par.chi_b+ ...
        (par.chi*eq.A*par.lambda)^(1/(1-par.alpha_q))*(par.alpha_q*eq.qw)^(par.alpha_q/(1-par.alpha_q))*...
        eq.lx^(par.alpha_x*par.alpha_q/(1-par.alpha_q)));
    f(2)=par.L_I-eq.lx-eq.lq;

end

%Numerically scale the equations
f=100*f;

end