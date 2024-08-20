%Solving the planners problem
%Solves for the labor allocation of the Social Planner's problem maximizing growth
function sp=sp_intq_B(par)

%Compute the optimal labor for quality
fun_lq=@(lq) (par.alpha_e*par.chi_e*par.lambda_e/(par.chi*par.lambda*par.alpha_q^(1-par.alpha_x)*par.alpha_x^par.alpha_x))^(1/(1-par.alpha_e))*lq.^((1-par.alpha_x-par.alpha_q)/(1-par.alpha_e))+...
           ((par.alpha_x+par.alpha_q)/par.alpha_q)*lq-par.L_I;
%lq0=0.1*par.L_I;
sp.lq=fzero(fun_lq,[0,par.L_I]);


%Labor for speed
sp.lx=(par.alpha_x/par.alpha_q)*sp.lq;

%Labor for entry
sp.le=par.L_I-sp.lx-sp.lq;

%Speed, quality and entry
sp.x=par.chi*sp.lx^par.alpha_x;
sp.Q=par.lambda*sp.lq^par.alpha_q;
sp.xe=par.chi_e*sp.le^par.alpha_e;

%Growth
sp.g=sp.x*sp.Q+sp.xe*par.lambda_e;

% %Check FOCs
% par.alpha_x*par.chi*sp.lx^(par.alpha_x-1)*par.lambda*sp.lq^par.alpha_q
% par.alpha_q*par.chi*sp.lx^par.alpha_x*par.lambda*sp.lq^(par.alpha_q-1)
% par.alpha_e*par.chi_e*sp.le^(par.alpha_e-1)*par.lambda_e
            
end