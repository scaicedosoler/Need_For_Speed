%Score function to estimate the parameters
function [f,par,eq,mom]=score_mom_intq_B(parv,dmom,par)

%Unpack parameters
if ~isfield(par,'par_est')
    par.par_est='par1'; %Make default
end

if strcmp(par.par_est,'par1')

par.lambda=parv(1);
par.chi_e=parv(2);
par.lambda_e=parv(3);
par.alpha_e=parv(4);
par.gamma_q=parv(5);
par.alpha_x=parv(7);
par.chi=parv(6);
par.chi_b=parv(7);

%To have simple solution of the model
par.gamma_e=1;

end

if strcmp(par.par_est,'par_nu')

par.lambda=parv(1);
par.chi_e=parv(2);
par.lambda_e=parv(3);
par.alpha_e=parv(4);
par.nu=parv(5);
par.gamma_q=parv(6);
par.alpha_x=parv(7);
par.chi=parv(8);
par.chi_b=parv(9);

end

if strcmp(par.par_est,'fix_chis')

par.lambda=parv(1);
par.chi_e=parv(2);
par.lambda_e=parv(3);
par.alpha_e=parv(4);
par.nu=parv(5);
par.gamma_q=parv(6);
par.alpha_x=parv(7);

end

if strcmp(par.par_est,'par_alpha_q')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.alpha_x=parv(6);
    par.alpha_q=parv(7);

end

if strcmp(par.par_est,'par_alpha_q_chi')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.alpha_x=parv(6);
    par.alpha_q=parv(7);
    par.chi=parv(8);

end

if strcmp(par.par_est,'par_alpha_q_chis')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.alpha_x=parv(6);
    par.alpha_q=parv(7);
    par.chi=parv(8);
    par.chi_b=parv(9);

end

if strcmp(par.par_est,'par_alpha_q_nu')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.alpha_x=parv(6);
    par.alpha_q=parv(7);
    par.nu=parv(8);

end

if strcmp(par.par_est,'par_all')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.alpha_x=parv(6);
    par.alpha_q=parv(7);
    par.chi=parv(8);
    par.chi_b=parv(9);
    par.nu=parv(10);

end

if strcmp(par.par_est,'par_chis_nu')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.chi=parv(6);
    par.chi_b=parv(7);
    par.nu=parv(8);

end

if strcmp(par.par_est,'par_no_alphas')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.chi=parv(6);
    par.chi_b=parv(7);
    par.nu=parv(8);

end

if strcmp(par.par_est,'par_L_I')

    par.lambda=parv(1);
    par.chi_e=parv(2);
    par.lambda_e=parv(3);
    par.alpha_e=parv(4);
    par.gamma_q=parv(5);
    par.chi=parv(6);
    par.chi_b=parv(7);
    par.nu=parv(8);
    par.L_I=parv(9);

end

if strcmp(par.par_est,'par_chib0')

%Initial parameters
par.lambda=parv(1);
par.chi_e=parv(2);
par.lambda_e=parv(3);
par.alpha_e=parv(4);
par.gamma_q=parv(5);
par.chi=parv(6);
par.nu=parv(7);
par.L_I=parv(8);


if par.iyear~=1980
    %Let chi_b change
    par.chi_b=parv(9);
end

end

if strcmp(par.par_est,'par_chib_50')

%Initial parameters
par.lambda=parv(1);
par.chi_e=parv(2);
par.lambda_e=parv(3);
par.alpha_e=parv(4);
par.gamma_q=parv(5);
par.chi=parv(6);
par.nu=parv(7);
par.L_I=parv(8);
par.chi_b=parv(9);

end

%Adjust parametric restrictions
par.gamma_x=1-par.alpha_x-par.alpha_q-par.gamma_q;
par.bargamma_x=par.alpha_x+par.alpha_q-par.bargamma_q;

par.gamma_b=par.gamma_q+par.alpha_q;
par.bargamma_b=par.bargamma_q+par.alpha_q;

par.xi_e=1-par.alpha_e-par.gamma_e;
par.barxi_e=par.alpha_e-par.bargamma_e;

%% Model moments

%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium
%%%%%%%%%%%%%%%%%%%%%%%

eq=eq_sim_fun_intq_B(par);

if isfield(eq, "err")
    f=100*ones(1,par.nmom); %Something large
    mom=[];
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%
% Moments
%%%%%%%%%%%%%%%%%%%%%%%

%1. Growth rate
mom.g=eq.g;

%2.Inventors in entrant firms
mom.inv_ent=eq.sle;

%3. Entrant innovations
mom.pat_ent=eq.se;

%4. Relative patent quality Q_e/Q
mom.q_ent_inc=eq.q_ent_inc;

%5. Inventors in top 10
mom.inv_10=eq.sl_top;

%6. Top 10 innovations
mom.pat_10=eq.s_top;

%7. Relative patent quality Q_top/Q_bot
mom.q_10_90=eq.q_top_bot;

%8. Patents relative to baseline
mom.pat_f_rb=eq.xbar/par.xbarb;

%9. Proportion of entry: top10 relative to bot90
mom.ent_10_90=(eq.xe_top/eq.x_top)/(eq.xe_bot/eq.x_bot);

%10. Average private value of innovations: (V/x top10)/(V/x bot 90) 
%mom.xi_10_90=(eq.V_top/eq.x_top)/(eq.V_bot/eq.x_bot);
mom.xi_10_90=eq.V_top/eq.V_bot;

%11. Ratio: private value/ patent quality relative to baseline
mom.rb_lxi_lq=(eq.Vbar/eq.dq)/par.V_dq_b;

%12. Patents per inventor
mom.x_l=(eq.x+eq.xe)/par.L_I;

%13. Ratio: private value/ patent quality relative to baseline
mom.rb_xi_lq=(eq.Vbar/eq.dq)/par.V_dq_b;

%14. Inventors per patent
mom.rb_inv_pat=(par.L_I/(eq.x+eq.xe))/par.inv_pat_b;

%15. Patent value over sale
mom.pat_val_sales=(eq.x*(eq.A*eq.dq+par.chi_b)/(eq.r+eq.tau))/eq.yj;

%16. Change in Patent value over sale relative to baseline
if ~isfield(par,'pat_val_sales_b')
    par.pat_val_sales_b=mom.pat_val_sales;
end

mom.rb_pat_val_sales=mom.pat_val_sales/par.pat_val_sales_b;

%17. Additional benefit over total profits
mom.pib_pi_tot=eq.pib/(eq.pi+eq.pib+eq.piq);

%% Score function

f=zeros(1,par.nmom);

f(1)=(mom.g-dmom.g)^2/(0.5*(mom.g)^2+0.5*dmom.g^2);
f(2)=(mom.inv_ent-dmom.inv_ent)^2/(0.5*(mom.inv_ent)^2+0.5*dmom.inv_ent^2);
f(3)=(mom.pat_ent-dmom.pat_ent)^2/(0.5*(mom.pat_ent)^2+0.5*dmom.pat_ent^2);
f(4)=(mom.q_ent_inc-dmom.q_ent_inc)^2/(0.5*(mom.q_ent_inc)^2+0.5*dmom.q_ent_inc^2);
f(5)=(mom.inv_10-dmom.inv_10)^2/(0.5*(mom.inv_10)^2+0.5*dmom.inv_10^2);
f(6)=(mom.pat_10-dmom.pat_10)^2/(0.5*(mom.pat_10)^2+0.5*dmom.pat_10^2);
f(7)=(mom.q_10_90-dmom.q_10_90)^2/(0.5*(mom.q_10_90)^2+0.5*dmom.q_10_90^2);
f(8)=(mom.pat_f_rb-dmom.pat_f_rb)^2/(0.5*(mom.pat_f_rb)^2+0.5*dmom.pat_f_rb^2);
f(9)=(mom.ent_10_90-dmom.ent_10_90)^2/(0.5*(mom.ent_10_90)^2+0.5*dmom.ent_10_90^2);
f(10)=(mom.xi_10_90-dmom.xi_10_90)^2/(0.5*(mom.xi_10_90)^2+0.5*dmom.xi_10_90^2);
f(11)=(mom.rb_lxi_lq-dmom.rb_lxi_lq)^2/(0.5*(mom.rb_lxi_lq)^2+0.5*dmom.rb_lxi_lq^2);
f(12)=(mom.x_l-dmom.x_l)^2/(0.5*(mom.x_l)^2+0.5*dmom.x_l^2);
f(13)=(mom.rb_xi_lq-dmom.rb_xi_lq)^2/(0.5*(mom.rb_xi_lq)^2+0.5*dmom.rb_xi_lq^2);
f(14)=(mom.rb_inv_pat-dmom.rb_inv_pat)^2/(0.5*(mom.rb_inv_pat)^2+0.5*dmom.rb_inv_pat^2);
f(15)=(mom.pat_val_sales-dmom.pat_val_sales)^2/(0.5*(mom.pat_val_sales)^2+0.5*dmom.pat_val_sales^2);
f(16)=(mom.rb_pat_val_sales-dmom.rb_pat_val_sales)^2/(0.5*(mom.rb_pat_val_sales)^2+0.5*dmom.rb_pat_val_sales^2);
f(17)=(mom.pib_pi_tot-dmom.pib_pi_tot)^2/(0.5*(mom.pib_pi_tot)^2+0.5*dmom.pib_pi_tot^2);

%Adjust units
f=sqrt(f);

end
