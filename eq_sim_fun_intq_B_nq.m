%Computes equilibrium of the model with internal innovation and additional
%benefit no endogenous quality
function [eq,sim]=eq_sim_fun_intq_B_nq(par)
%% Parameters and options

%Specific parameter restrictions
par.gamma_x=1-par.alpha_x-par.alpha_q-par.gamma_q;
par.bargamma_x=par.alpha_x+par.alpha_q-par.bargamma_q;

par.gamma_b=par.gamma_q+par.alpha_q;
par.bargamma_b=par.bargamma_q+par.alpha_q;

par.xi_e=1-par.alpha_e-par.gamma_e;
par.barxi_e=par.alpha_e-par.bargamma_e;

%Running options
sol_opt='fm_con'; %'fs'; % 'fm';
run_sim=0;

%Plot distribution
plot_dist=0;

%If no parameter options for algorithm

if isempty(par.alg.opt_fs)
    %Options for solvers
    par.alg.opt_fz = optimset('Display','off');
    par.alg.opt_fs = optimoptions('fsolve','Display','off');
    par.alg.opt_fmincon0=optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
    par.alg.opt_fmincon=optimoptions(@fmincon,'Algorithm','interior-point','Display','off','OutputFcn',@stopfun);
end


%% Solving the equilibrium

fun_var=@(var)sys_eq_sim_intq_B_nq(var,par);

var0=0.6*par.L_I;

fun_var(var0);

%Lower and upper bounds
lb = 0;
ub= par.L_I;

%Linear inequality constraints
A=[];
b=[];

%Equality constraints
Ae=[];
be=[];

%Tolerance
tol=1e-4;


if strcmp(sol_opt,'fs')

    %----------------------
    % Solve using fsolve
    %----------------------

    tic
    [var, fval, exitflag]=fsolve(fun_var,var0,par.alg.opt_fs);
    toc

    if exitflag>0 && sum(abs(fval))<tol
        disp('fsolve converged')

        %Equilibrium variables
            [res,eq]=sys_eq_sim_intq_B_nq(var,par);

        eq.res=res;
        eq.exitflag=exitflag;

        return

    else
        disp('fsolve: No convergence')
        sol_opt='fm';
    end

end

if strcmp(sol_opt,'fm')
    %----------------------
    % Solve using fmincon
    %----------------------
    disp('Using fmincon to minimize sum of squared errors...')

    error_fun=@(var) sum(fun_var(var).^2); %Sum of squared errors


    %--------------------------------
    % Solve using fmincon
    %--------------------------------
    tic
    [var,fval,exitflag] = fmincon(error_fun,var0,A,b,Ae,be,lb,ub,[],par.alg.opt_fmincon0);
    toc

    if exitflag>0 && abs(fval)<tol
        disp('fmincon converged')

        %Equilibrium variables
            [res,eq]=sys_eq_sim_intq_B_nq(var,par);

        eq.res=res;
        eq.exitflag=exitflag;

        %return
    else
        disp('fmincon: No convergence')
        sol_opt='fm_con';
    end

end

if strcmp(sol_opt,'fm_con')
    %-------------------------------------
    % Solve using fmincon as constraints
    %-------------------------------------

    disp('Using fmincon to solve equilibrium equations as constraints...')

    tic
    [var,fval,exitflag] = fmincon(@(var)0,var0,A,b,Ae,be,lb,ub,@(var)fminconstr_intq_B_nq(var,par),par.alg.opt_fmincon);
    toc

    %Equilibrium variables
    [res,eq]=sys_eq_sim_intq_B_nq(var,par);


    eq.res=res;

    if exitflag>0 && isfield(eq,"x")
        disp('fmincon constraints converged')
        eq.exitflag=exitflag;

        %return
    else
        disp('fmincon constraints: No convergence')
        eq=[];
        sim=[];
        eq.res=res;
        eq.exitflag=exitflag;

        eq.err=1;

        return
    end

end


%% Simulations : Should update--Mar 2024
% 
% if run_sim==1
%     %% Run simulations
% 
%     %To test
%     %eq.xe=0;
%     %par.alpha_q=1;
%     disp('Running simulations...')
% 
%     tic
%         sim=sim_fun_intq_B_nu(par,eq);
%     toc
% 
% 
%     %Distribution
%     eq.Tsim=sim.Tsim;
%     eq.qsim=sim.Q(:,sim.Tsim);
%     [sim.FTsim,sim.QTsim] = ecdf(sim.Q(:,sim.Tsim));
% 
%     eq.QTsim=sim.QTsim;
%     eq.FTsim=sim.FTsim;
% 
%     %Divide incumbents into top10 and bot90
%     sim.Q_p90=prctile(sim.Q(:,sim.Tsim),90);
%     sim.IQtop10=(sim.Q(:,sim.Tsim)>sim.Q_p90);
%     sim.IQbot90=~sim.IQtop10;
% 
%     % %Tail of the distribution: Log(1-F(q))
%     % etail=1e-2;
%     % ftail=1-etail;
%     % Itail=find(sim.FTsim>=ftail&sim.FTsim<1);
%     %
%     % sim.Ftailq=log(1-sim.FTsim(Itail));
%     % sim.logqtail=log(sim.QTsim(Itail));
%     %
%     % %Fit a regression
%     % eq.coef = [ones(length(sim.logqtail),1) sim.logqtail]\sim.Ftailq;
%     % eq.tailcons=sim.Ftailq(1)-eq.coef(2)*sim.logqtail(1);
%     % eq.tailhat=eq.tailcons+eq.coef(2)*sim.logqtail;
%     %
%     % eq.tailpar=-1/eq.coef(2);
% 
%     %Compute the average quality for top10 and bot90
%     eq.Qtop=mean(sim.Q(sim.IQtop10,sim.Tsim));
%     eq.Qbot=mean(sim.Q(sim.IQbot90,sim.Tsim));
%     eq.Qbar=mean(sim.Q(:,sim.Tsim));
% 
%     if plot_dist==1
% 
%         %Quality distribution-product line
%         tp25=floor((sim.Tsim-1)/4);
%         tp50=floor((sim.Tsim-1)/2);
%         tp75=floor(3*(sim.Tsim-1)/4);
% 
%         [sim.F0,sim.Q0] = ecdf(sim.Q(:,1));
%         [sim.Fp25,sim.Qp25] = ecdf(sim.Q(:,tp25));
%         [sim.Fp50,sim.Qp50] = ecdf(sim.Q(:,tp50));
%         [sim.Fp75,sim.Qp75] = ecdf(sim.Q(:,tp75));
% 
%         % Plot growth
%         figure; hold all;
%         plot(sim.g/par.sim.dt);
%         plot(mean(sim.g/par.sim.dt)*ones(1,sim.Tsim),'--')
%         plot(eq.g*ones(1,sim.Tsim),'-k')
%         xlabel('Period')
%         ylabel('Growth')
% 
%         %Distribution of product lines quality
%         figure; hold all;
%         plot([0,0], [0,1],':k');
%         p0=plot(sim.Q0,sim.F0,'--k');
%         p25=plot(sim.Qp25,sim.Fp25,'--');
%         p50=plot(sim.Qp50,sim.Fp50,'--');
%         p75=plot(sim.Qp75,sim.Fp75,'--');
%         pT=plot(sim.QTsim,sim.FTsim); pT.Color=par.opt.light_blue;
%         lgd=legend([p0,p25,p50,p75,pT],{'Initial','$t_{p25}$','$t_{p50}$','$t_{p75}$', 'Final'},'Interpreter','latex'); lgd.Location='best'; lgd.EdgeColor='none';
%         xlim([0,10])
%         xlabel('Quality')
%         ylabel('CDF')
% 
%         % Concentration
%         figure; hold all;
%         plot(0.1*sim.QbarT./sim.Qbar);
%         ylim([0,1])
%         xlabel('Period')
%         ylabel('Concentration of Quality in Top 10')
% 
% 
%     end
% 
% 
%     %% Compute moments by groups
% 
%     %Compute the arrival rate distribution
%     sim.xq=par.chi*eq.lx^par.alpha_x*eq.qsim.^(par.alpha_x+par.gamma_x)*mean(eq.qsim).^(par.bargamma_x-par.alpha_x);
%     sim.Qq=par.lambda*eq.lq^par.alpha_q*eq.qsim.^(par.alpha_q+par.gamma_q).*mean(eq.qsim).^(par.bargamma_q-par.alpha_q);
% 
%     %Numerical adjustment of wages
%     eq.w=eq.w*(eq.Qbar/par.qbar);
% 
%     eq.qw=par.qbar/eq.w;
% 
%     %Numerical adjustment of Clx and Clq 
%     %eq=adj_Clx_Clq(par,eq);
% 
%     %Aggregate labor: speed and quality
%     eq.lx_top=(1-par.varrho)*eq.Phi*eq.Clx*(eq.Qtop/eq.Qbar);
%     eq.lq_top=(1-par.varrho)*eq.Phi*eq.Clq*(eq.Qtop/eq.Qbar);
% 
%     eq.lx_bot=par.varrho*eq.Phi*eq.Clx*(eq.Qbot/eq.Qbar);
%     eq.lq_bot=par.varrho*eq.Phi*eq.Clq*(eq.Qbot/eq.Qbar);
% 
% 
%     % Arrival rate
%     eq.qxtop=mean(sim.Q(sim.IQtop10,sim.Tsim).^(par.alpha_x+par.gamma_x));
%     eq.qxbot=mean(sim.Q(sim.IQbot90,sim.Tsim).^(par.alpha_x+par.gamma_x));
%     eq.qxbar=mean(sim.Q(:,sim.Tsim).^(par.alpha_x+par.gamma_x));
% 
%     eq.x_top=(1-par.varrho)*par.chi*eq.Clx^par.alpha_x*eq.Qbar^(par.bargamma_x-par.alpha_x)*eq.qxtop;
%     eq.x_bot=par.varrho*par.chi*eq.Clx^par.alpha_x*eq.Qbar^(par.bargamma_x-par.alpha_x)*eq.qxbot;
%     eq.xbar=par.chi*eq.Clx^par.alpha_x*eq.Qbar^(par.bargamma_x-par.alpha_x)*eq.qxbar;
% 
% 
%     sim.Qq=par.lambda*eq.lq^par.alpha_q*sim.Q(:,sim.Tsim).^(par.alpha_q+par.gamma_q).*mean(sim.Q(:,sim.Tsim)).^(par.bargamma_q-par.alpha_q);
% 
%     %Quality
%     eq.qQtop=mean(sim.Q(sim.IQtop10,sim.Tsim).^(par.alpha_q+par.gamma_q));
%     eq.qQbot=mean(sim.Q(sim.IQbot90,sim.Tsim).^(par.alpha_q+par.gamma_q));
%     eq.qQbar=mean(sim.Q(:,sim.Tsim).^(par.alpha_q+par.gamma_q));
% 
%     eq.Q_top=par.lambda*eq.Clq^par.alpha_q*eq.Qbar^(par.bargamma_q-par.alpha_q)*eq.qQtop;
%     eq.Q_bot=par.lambda*eq.Clq^par.alpha_q*eq.Qbar^(par.bargamma_q-par.alpha_q)*eq.qQbot;
%     eq.Q_bar=par.lambda*eq.Clq^par.alpha_q*eq.Qbar^(par.bargamma_q-par.alpha_q)*eq.qQbar;
% 
%     %Additional benefit
%     eq.qBtop=mean(sim.Q(sim.IQtop10,sim.Tsim).^par.gamma_b);
%     eq.qBbot=mean(sim.Q(sim.IQbot90,sim.Tsim).^par.gamma_b);
%     eq.qBbar=mean(sim.Q(:,sim.Tsim).^par.gamma_b);
% 
%     eq.B_top=(1-par.varrho)*par.chi_b*eq.Qbar^par.bargamma_b*eq.qBtop;
%     eq.B_bot=par.varrho*par.chi_b*eq.Qbar^par.bargamma_b*eq.qBbot;
%     eq.B_bar=par.chi_b*eq.Qbar^par.bargamma_b*eq.qBbar;
% 
%     %Entry
%     eq.Cle=(par.alpha_e*par.chi_e*eq.A*(1+par.lambda_e))^(1/(1-par.alpha_e));
% 
%     % eq.qxetop=mean(sim.Q(sim.IQtop10,sim.Tsim).^(par.alpha_e+par.xi_e));
%     % eq.qxebot=mean(sim.Q(sim.IQbot90,sim.Tsim).^(par.alpha_e+par.xi_e));
%     % eq.qxebar=mean(sim.Q(:,sim.Tsim).^(par.alpha_e+par.xi_e));
%     % 
%     % eq.xe_top=(1-par.varrho)*par.chi_e*eq.Cle^par.alpha_e*eq.Qbar^(par.barxi_e-par.alpha_e)*eq.qxetop;
%     % eq.xe_bot=par.varrho*par.chi_e*eq.Cle^par.alpha_e*eq.Qbar^(par.barxi_e-par.alpha_e)*eq.qxebot;
%     % eq.xebar=par.chi_e*eq.Cle^par.alpha_e*eq.Qbar^(par.barxi_e-par.alpha_e)*eq.qxbar;
% 
%     %Check this!
%     eq.qxetop=mean(sim.Q(sim.IQtop10,sim.Tsim));
%     eq.qxebot=mean(sim.Q(sim.IQbot90,sim.Tsim));
%     eq.qxebar=mean(sim.Q(:,sim.Tsim));
% 
%     eq.xe_top=(1-par.varrho)*par.chi_e*eq.Cle^par.alpha_e*eq.Qbar*eq.qxetop;
%     eq.xe_bot=par.varrho*par.chi_e*eq.Cle^par.alpha_e*eq.Qbar*eq.qxebot;
%     eq.xebar=par.chi_e*eq.Cle^par.alpha_e*eq.Qbar*eq.qxbar;
% 
%     % Aggregate labor top 10 and bottom 90 firms
%     eq.l_top=eq.lx_top+eq.lq_top;
%     eq.l_bot=eq.lx_bot+eq.lq_bot;
% 
%     %Private value
%     eq.V_top=(1-par.varrho)*eq.A*eq.Qtop;
%     eq.V_bot=par.varrho*eq.A*eq.Qbot;
%     eq.Vbar=eq.A*eq.Qbar;
% 
% 
% 
%     %% Moments
% 
%     comp_mom=1;
% 
%     if comp_mom==1
% 
%         %Growth
%         %eq.g=eq.g;
% 
%         %--------------------
%         % Share of inventors
%         %--------------------
% 
%         % Inventors in top 10
%         eq.sl_top=eq.l_top/par.L_I;
% 
%         % Inventors in entrant firms
%         eq.sle=eq.le/par.L_I;
% 
%         %----------------------
%         % Share of innovations
%         %----------------------
% 
%         % Innovations in top 10
%         eq.s_top=eq.x_top/(eq.xbar+eq.xe);
% 
%         % Entrant innovations
%         eq.se=eq.xe/(eq.x+eq.xe);
% 
% 
%         %----------------------
%         % Patent quality
%         %----------------------
%         eq.q_ent_inc=eq.dqe/eq.Q_bar;
%         eq.q_top_bot=eq.Q_top/eq.Q_bot;
% 
% 
% 
%         %%%%%%%%%%%%%%%%%%%%%%%
%         % Moments (in structure)
%         %%%%%%%%%%%%%%%%%%%%%%%
% 
%         %1. Growth rate
%         eq.mom.g=eq.g;
% 
%         %2.Inventors in entrant firms
%         eq.mom.inv_ent=eq.sle;
% 
%         %3. Entrant innovations
%         eq.mom.pat_ent=eq.se;
% 
%         %4. Relative patent quality Q_e/Q
%         eq.mom.q_ent_inc=eq.q_ent_inc;
% 
%         %5. Inventors in top 10
%         eq.mom.inv_10=eq.sl_top;
% 
%         %6. Top 10 innovations
%         eq.mom.pat_10=eq.s_top;
% 
%         %7. Relative patent quality Q_top/Q_bot
%         eq.mom.q_10_90=eq.q_top_bot;
% 
%         %8. Patents relative to baseline
%         eq.mom.pat_rb=eq.xbar/par.xbarb;
% 
%         %9. Proportion of entry: top10 relative to bot90  
%         eq.mom.ent_10_90=(eq.xe_top/eq.x_top)/(eq.xe_bot/eq.x_bot);
% 
%         %10. Average private value of innovations: (V/x top10)/(V/x bot 90) 
%         %eq.mom.xi_10_90=(eq.V_top/eq.x_top)/(eq.V_bot/eq.x_bot);
%         eq.mom.xi_10_90=eq.V_top/eq.V_bot;
% 
%         %11. Ratio: private value/ patent quality relative to baseline
%         eq.mom.r_lxi_lq=(eq.Vbar/eq.dq)/par.V_dq_b;
% 
%         %12. Patents per inventor
%         eq.mom.x_l=(eq.x+eq.xe)/par.L_I;
% 
% 
%         %--------------------------------------------
%         %--------------------------------------------
%         % Extended Moments
%         %--------------------------------------------
%         %--------------------------------------------
% 
%         comp_ext_mom=1;
% 
%         if comp_ext_mom==1
% 
%         %Functions
%         [qsim_u, ind_u]=unique(eq.qsim); %Make sure the inperpolation is well-defined
% 
%         x_fun=@(q)interp1(qsim_u,sim.xq(ind_u),q,"linear","extrap");
%         dq_fun=@(q)interp1(qsim_u,sim.Qq(ind_u),q,"linear","extrap");
% 
%         % Firm Size by percentiles
%         %Divide incumbents into percentile ranges
%         eq.Q_p50=prctile(eq.qsim,50);
%         eq.Q_p75=prctile(eq.qsim,75);
%         eq.Q_p90=prctile(eq.qsim,90);
%         eq.Q_p99=prctile(eq.qsim,99);
%         eq.Q_p999=prctile(eq.qsim,99.9);
% 
%         sim.IQp50=(eq.qsim<=eq.Q_p50);
%         sim.IQp75_50=(eq.Q_p50<eq.qsim&eq.qsim<=eq.Q_p75);
%         sim.IQp90_75=(eq.Q_p75<eq.qsim&eq.qsim<=eq.Q_p90);
%         sim.IQp99_90=(eq.Q_p90<eq.qsim&eq.qsim<eq.Q_p99);
%         sim.IQp99a=(eq.Q_p99<=eq.qsim);
% 
%         %Adjust top percentiles
%         if mean(sim.IQp99a)>0.01
% 
%             %Identify top 1%
%             idx=find(sim.IQp99a==1);
%             idx99a=randsample(idx,100);
%             idx90_99= setdiff(idx,idx99a);
%             %indx99=indx(1:100); %Faster but less general
% 
%             sim.IQp99a(idx)=0;
%             sim.IQp99a(idx99a)=1;
% 
%             sim.IQp99_90(idx90_99)=1;
% 
%         end
% 
%         %Average quality level for each level
%         eq.Qp50bar=mean(sim.Q(sim.IQp50,sim.Tsim));
%         eq.Qp75_50bar=mean(sim.Q(sim.IQp75_50,sim.Tsim));
%         eq.Qp90_75bar=mean(sim.Q(sim.IQp90_75,sim.Tsim));
%         eq.Qp99_90bar=mean(sim.Q(sim.IQp99_90,sim.Tsim));
%         eq.Qp99abar=mean(sim.Q(sim.IQp99a,sim.Tsim));
% 
%         %-------------------
%         %Moments
%         %-------------------
%         %Weights of percentile groups
%         eq.pg={'p50', '(p50,p75]', '(p75,p90]', '(p90,p99]', 'Above p99'};
%         eq.npg=length(eq.pg);
%         eq.vpg=1:eq.npg;
%         eq.wpg=[0.5,0.25,0.15,0.09,0.01];
% 
%         %Speed and entry
% 
%         %Average arrival rate by percentile group
%         eq.x_pg=[mean(max(x_fun(sim.Q(sim.IQp50,sim.Tsim)),0)),...
%             mean(max(x_fun(sim.Q(sim.IQp75_50,sim.Tsim)),0)),...
%             mean(max(x_fun(sim.Q(sim.IQp90_75,sim.Tsim)),0)),...
%             mean(max(x_fun(sim.Q(sim.IQp99_90,sim.Tsim)),0)),...
%             mean(max(x_fun(sim.Q(sim.IQp99a,sim.Tsim)),0))];
% 
%         eq.x_wpg=eq.x_pg.*eq.wpg;
% 
%         eq.xe_pg=[eq.xe*eq.Qp50bar, eq.xe*eq.Qp75_50bar, eq.xe*eq.Qp90_75bar, eq.xe*eq.Qp99_90bar, eq.xe*eq.Qp99abar];
% 
%         eq.xe_wpg=eq.xe_pg.*eq.wpg;
% 
% 
%         eq.lnx_pg=[mean(log(max(x_fun(sim.Q(sim.IQp50,sim.Tsim)),eps))),...
%             mean(log(max(x_fun(sim.Q(sim.IQp75_50,sim.Tsim)),eps))),...
%             mean(log(max(x_fun(sim.Q(sim.IQp90_75,sim.Tsim)),eps))),...
%             mean(log(max(x_fun(sim.Q(sim.IQp99_90,sim.Tsim)),eps))),...
%             mean(log(max(x_fun(sim.Q(sim.IQp99a,sim.Tsim)),eps)))];
% 
%         eq.lnx_wpg=eq.lnx_pg.*eq.wpg;
% 
%         eq.lnxe_pg=log(eq.xe_pg);
% 
%         eq.lnxe_wpg=eq.lnxe_pg.*eq.wpg;
% 
% 
%         %Labor
%         eq.lx_pg=[eq.lx*eq.Qp50bar, eq.lx*eq.Qp75_50bar, eq.lx*eq.Qp90_75bar, eq.lx*eq.Qp99_90bar, eq.lx*eq.Qp99abar];
% 
%         eq.lx_wpg=eq.lx_pg.*eq.wpg;
% 
%         eq.lq_pg=[eq.lq*eq.Qp50bar, eq.lq*eq.Qp75_50bar, eq.lq*eq.Qp90_75bar, eq.lq*eq.Qp99_90bar, eq.lq*eq.Qp99abar];
% 
%         eq.lq_wpg=eq.lq_pg.*eq.wpg;
% 
%         eq.le_pg=[eq.le*eq.Qp50bar, eq.le*eq.Qp75_50bar, eq.le*eq.Qp90_75bar, eq.le*eq.Qp99_90bar, eq.le*eq.Qp99abar];
% 
%         eq.le_wpg=eq.le_pg.*eq.wpg;
% 
%         %Pantent quality
%         eq.dq_pg=[mean(max(dq_fun(sim.Q(sim.IQp50,sim.Tsim)),0)),...
%             mean(max(dq_fun(sim.Q(sim.IQp75_50,sim.Tsim)),0)),...
%             mean(max(dq_fun(sim.Q(sim.IQp90_75,sim.Tsim)),0)),...
%             mean(max(dq_fun(sim.Q(sim.IQp99_90,sim.Tsim)),0)),...
%             mean(max(dq_fun(sim.Q(sim.IQp99a,sim.Tsim)),0))];
% 
%         %eq.dqe_pg=[];
% 
%         %Firm Value
%         eq.V_pg=[eq.A*eq.Qp50bar, eq.A*eq.Qp75_50bar, eq.A*eq.Qp90_75bar, eq.A*eq.Qp99_90bar, eq.A*eq.Qp99abar];
%         eq.Vx_pg=eq.V_pg./eq.x_pg;
% 
%         end
%     end
% 
% end




end