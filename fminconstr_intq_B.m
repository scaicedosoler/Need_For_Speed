function [c,ceq] = fminconstr_intq_B(var,par)

c = []; % No nonlinear inequality
ceq = sys_eq_sim_intq_B(var,par); %  Nonlinear equality constraints

end