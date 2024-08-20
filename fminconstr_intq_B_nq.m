function [c,ceq] = fminconstr_intq_B_nq(var,par)

c = []; % No nonlinear inequality

% Solve system of equations for var=(r,w,le)
ceq = sys_eq_sim_intq_B_nq(var,par); %  Nonlinear equality constraints


end