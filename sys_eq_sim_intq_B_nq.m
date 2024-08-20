%System of equations to find equilibrium with internal innovation and additional
%benefit labor--no endogenous quality
function [f,eq]=sys_eq_sim_intq_B_nq(var,par)

%Unpack variables var=(lx)
eq.lx=var(1);

%Mass of product lines
eq.Phi=1; 

% Minimum quality
eq.qmin=0; 

%Benefit
eq.B=par.chi_b;

%Production
eq.betatilde=par.beta^par.beta*(1-par.beta)^(1-2*par.beta);
eq.lp=par.beta*par.Lp/(par.beta+eq.Phi*(1-par.beta)^2);
eq.pi=eq.lp*(1-par.beta)*eq.betatilde;

eq.Y=(eq.betatilde/par.beta)*par.qbar*eq.lp;
eq.wp=eq.betatilde*par.qbar;
eq.pj=eq.wp/((1-par.beta)*par.qbar);
eq.kj=eq.lp*((1-par.beta)*par.qbar/eq.wp).^(1/par.beta);
eq.yj=eq.kj*eq.pj;


%Arrival rate
eq.x=par.chi*eq.lx^par.alpha_x;

%Labor 
eq.Clx=eq.lx;

%Quality
eq.Q=par.lambda;

%Entrants
eq.le=par.L_I-eq.lx;

%Creative destruction
eq.xe=par.chi_e*eq.le^par.alpha_e;
eq.tau=eq.xe;

%Quality
eq.dq=par.lambda;
eq.dqe=par.lambda_e;

%Growth
eq.g=eq.xe*eq.dqe+eq.x*eq.dq;

%Interest rate
eq.r=eq.g+par.rho;

%8.A 
eq.pib=par.chi*par.chi_b*eq.Clx^(par.alpha_x)*(1-par.alpha_x);
eq.piq=0;

eq.A=(eq.pi+eq.pib)/((eq.r+eq.tau)-par.chi*par.lambda*eq.Clx^(par.alpha_x)*(1-par.alpha_x));

%9.w
eq.qw=eq.le^(1-par.alpha_e)/(par.alpha_e*par.chi_e*eq.A*(1+par.lambda_e));
eq.wq=1/eq.qw;
eq.w=eq.wq/par.qbar;

%--------------------------
% Equations
%--------------------------

f(1)=eq.lx^(1-par.alpha_x)-par.alpha_x*eq.qw*(par.chi*par.chi_b+ ...
    (par.chi*eq.A*par.lambda)^(1/(1-par.alpha_q))*(par.alpha_q*eq.qw)^(par.alpha_q/(1-par.alpha_q))*eq.lx^(par.alpha_x*par.alpha_q/(1-par.alpha_q)));

%Numerically scale the equations
f=100*f;

end