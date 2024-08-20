%Score function to estimate the parameters
function [f,par,eq,mom]=score_mom_intq_B_nq(parv,dmom,par)

%Unpack parameters
par.lambda=parv(1);

% %Adjust parametric restrictions
% par.gamma_x=1-par.alpha_x-par.alpha_q-par.gamma_q;
% par.bargamma_x=par.alpha_x+par.alpha_q-par.bargamma_q;
% 
% par.gamma_b=par.gamma_q+par.alpha_q;
% par.bargamma_b=par.bargamma_q+par.alpha_q;
% 
% par.xi_e=1-par.alpha_e-par.gamma_e;
% par.barxi_e=par.alpha_e-par.bargamma_e;

%% Model moments

%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium
%%%%%%%%%%%%%%%%%%%%%%%

eq=eq_sim_fun_intq_B_nq(par);

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

% %2.Inventors in entrant firms
% mom.inv_ent=eq.sle;
% 
% %3. Entrant innovations
% mom.pat_ent=eq.se;
% 
% %4. Relative patent quality Q_e/Q
% mom.q_ent_inc=eq.q_ent_inc;
%
% %5. Inventors in top 10
% mom.inv_10=eq.sl_top;
% 
% %6. Top 10 innovations
% mom.pat_10=eq.s_top;
% 
% %7. Relative patent quality Q_top/Q_bot
% mom.q_10_90=eq.q_top_bot;
% 
% %8. Patents relative to baseline
% mom.pat_rb=eq.xbar/par.xbarb;
% 
% %9. Proportion of entry: top10 relative to bot90
% mom.ent_10_90=(eq.xe_top/eq.x_top)/(eq.xe_bot/eq.x_bot);
% 
% %10. Average private value of innovations: (V/x top10)/(V/x bot 90) 
% %mom.xi_10_90=(eq.V_top/eq.x_top)/(eq.V_bot/eq.x_bot);
% mom.xi_10_90=eq.V_top/eq.V_bot;
% 
% %11. Ratio: private value/ patent quality relative to baseline
% mom.r_lxi_lq=(eq.Vbar/eq.dq)/par.V_dq_b;
% 
% %12. Patents per inventor
% mom.x_l=(eq.x+eq.xe)/par.L_I;
% 
% %% Adjust data moments
% 
% %Account for effect of growth on number of patents
% % g_adj=(1+(par.gamma_x+par.bargamma_x)*eq.g)^30;
% % dmom.pat_rb=dmom.pat_rb/g_adj;
% %Patents per firm
% dmom.pat_rb=dmom.pat_f_rb;
% 
% %Additional moments
% if par.iyear<=1995
%         %Arrival date over unique inventors
%         dmom.x_l=0.76; %Approximation from chart...must compute it correctly : Feb 2024
% 
% end
% 
% if par.iyear>1995
%         %Arrival date over unique inventors
%         dmom.x_l=0.67; %Approximation from chart...must compute it correctly : Feb 2024
% end


%% Score function

f=zeros(1,par.nmom);

f(1)=(mom.g-dmom.g)^2/(0.5*(mom.g)^2+0.5*dmom.g^2);
% f(2)=(mom.inv_ent-dmom.inv_ent)^2/(0.5*(mom.inv_ent)^2+0.5*dmom.inv_ent^2);
% f(3)=(mom.pat_ent-dmom.pat_ent)^2/(0.5*(mom.pat_ent)^2+0.5*dmom.pat_ent^2);
% f(4)=(mom.q_ent_inc-dmom.q_ent_inc)^2/(0.5*(mom.q_ent_inc)^2+0.5*dmom.q_ent_inc^2);
% f(5)=(mom.inv_10-dmom.inv_10)^2/(0.5*(mom.inv_10)^2+0.5*dmom.inv_10^2);
% f(6)=(mom.pat_10-dmom.pat_10)^2/(0.5*(mom.pat_10)^2+0.5*dmom.pat_10^2);
% f(7)=(mom.q_10_90-dmom.q_10_90)^2/(0.5*(mom.q_10_90)^2+0.5*dmom.q_10_90^2);
% f(8)=(mom.pat_rb-dmom.pat_rb)^2/(0.5*(mom.pat_rb)^2+0.5*dmom.pat_rb^2);
% f(9)=(mom.ent_10_90-dmom.ent_10_90)^2/(0.5*(mom.ent_10_90)^2+0.5*dmom.ent_10_90^2);
% f(10)=(mom.xi_10_90-dmom.xi_10_90)^2/(0.5*(mom.xi_10_90)^2+0.5*dmom.xi_10_90^2);
% f(11)=(mom.r_lxi_lq-dmom.r_lxi_lq)^2/(0.5*(mom.r_lxi_lq)^2+0.5*dmom.r_lxi_lq^2);
% f(12)=(mom.x_l-dmom.x_l)^2/(0.5*(mom.x_l)^2+0.5*dmom.x_l^2);

%Adjust units
f=sqrt(f);

end
