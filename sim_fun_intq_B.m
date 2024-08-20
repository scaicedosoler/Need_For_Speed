%Simulates model and computes average quality

function sim=sim_fun_intq_B(par,eq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seed=par.sim.seed;
nprod=par.sim.nprod;
Tmax=par.sim.Tmax;
Tmin=2000;
dt=par.sim.dt;

if ~isfield(par,'nu')
    par.nu=1;
end

% Set seed
rng(seed);

% Options
load_dist0=1; %load a previous quality and product line property realization
save_Qdist=0; %save distribution
max_sim=1; %use the maximum for quality improvement
random_sim=0; %take a random quality improvement


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Product line quality
sim.Q= zeros(nprod,Tmax); 
eq.qbar=1;

%Firm owning the product line
sim.FP= zeros(nprod,Tmax);

if load_dist0==1
    load('sim_intq_B_Q0.mat','Q0');
    sim.Q(:,1)=Q0;

else
    dist0='ll'; %''; %'exp'; % 
    F0=rand(nprod,1);

    %Exponential
    if strcmp(dist0,'exp')
        inv_F0_fun=@(F) -eq.qbar*log(1-F); %Suppose the initial distribution is exponential with mean = qbar
        sim.Q(:,1)=inv_F0_fun(F0);

    elseif strcmp(dist0,'ll')
        %Log-logistic
        theta0=0.5; %Thickness of the tail, theta=1 is Gibrat's law
        x0=(sin(pi*theta0)/(pi*theta0))^(1/theta0);
        inv_F0_fun=@(F) (x0*F./(1-F)).^theta0; %Suppose the initial distribution is log-logistic with mean = qbar
        sim.Q(:,1)=inv_F0_fun(F0);

    else
        %Uninformative
        sim.Q(:,1)=eq.qbar*ones(nprod,1);

    end
end

%Indicator of active product lines (above minimum quality) and of measure
sim.Iact=true(nprod,Tmax);
eq.qmin=0;
sim.Iact(:,1)=(sim.Q(:,1)>eq.qmin);

%Active lines are of measure eq.Phi
%sim.AL= binornd(1,eq.Phi,[nprod,Tmax]).*sim.Iact;
sim.AL= sim.Iact;

%Probability of successful innovation
sim.Px= zeros(nprod,Tmax);
sim.Pe= zeros(nprod,Tmax);

%External and entrant innovations
sim.xe_a(:,1)=min(eq.xe*dt,1); %Adjusted so probability is bounded
sim.Pe(:,1)= binornd(1,sim.xe_a(:,1),[nprod,1]);

%Internal innovations
sim.x_a(:,1)=min(eq.x*dt*sim.Q(:,1).^(par.alpha_x+par.gamma_x)*mean(sim.Q(sim.Iact(:,1),1)).^(par.bargamma_x-par.alpha_x),1); %Adjusted so probability is bounded
sim.Px(:,1)= binornd(1,sim.x_a(:,1),[nprod,1]);

%Run simulations
sim.DQ= zeros(nprod,Tmax);
sim.I= ones(nprod,Tmax);

%Compute distance of cdf percentiles
nprc=100;
prc=linspace(1,99,nprc);

sim.Q_prc=zeros(nprc,Tmax);
sim.Prc_Q=zeros(nprc,Tmax);

sim.Q_prc(:,1)=prctile(sim.Q(:,1),prc);
sim.Prc_Q(:,1)=invprctile(sim.Q(:,1),sim.Q_prc(:,1));


%Compute the average Q
sim.Qbar=zeros(1,Tmax);
sim.Qbar(1)=mean(sim.Q(sim.Iact(:,1),1));

sim.QbarT=zeros(1,Tmax);
sim.QbarB=zeros(1,Tmax);

%Residual of distribution
sim.res=ones(1,Tmax);
dt_res=20/dt; %Number of periods for rolling window to compute residual

sim.res_conc=1; sim.res_conp=0;
dt_con=100/dt; % Number of periods to run regression

%Iterations
it=1;

while (par.sim.tol<sim.res_conc && sim.res_conp<=0.05 && it<=Tmax)||it<=Tmin %The criteria using concentration

    %Update counter
    it=it+1;

    %Maximum Quality improvement
    if max_sim==1

    A=sim.Px(:,it-1).*par.lambda.*eq.lq.^(par.alpha_q).*(sim.Q(:,it-1)/sim.Qbar(it-1)).^(par.alpha_q+par.gamma_q).*sim.Qbar(it-1).*sim.AL(:,it-1);
    B=sim.Pe(:,it-1).*(par.lambda_e.*(par.nu*sim.Q(:,it-1)+(1-par.nu)*sim.Qbar(it-1)));
    sim.DQ(:,it)=double(A>B).*A+double(A<B).*B; %We get the maximum for each product line
    sim.I(:,it)=double(A>B)*1+double(A<B)*2+1; %We get the index for each product line

    end
    
    %Random Quality improvement
    if random_sim==1
    sim.I=randi([0, 1], nprod, it);
    A=sim.Px(:,it-1).*par.lambda.*eq.lq.^(par.alpha_q).*(sim.Q(:,it-1)/sim.Qbar(it-1)).^(par.alpha_q+par.gamma_q).*sim.Qbar(it-1).*sim.AL(:,it-1);
    B=sim.Pe(:,it-1).*(par.lambda_e.*(par.nu*sim.Q(:,it-1)+(1-par.nu)*sim.Qbar(it-1)));
    sim.DQ(:,it)=(A+B).*(1-sim.Px(:,it-1).*sim.Pe(:,it-1))+...
        sim.Px(:,it-1).*sim.Pe(:,it-1).*(sim.I(:,it).*A+(1-sim.I(:,it)).*B);

    end
       
    %Update quality
    sim.Q(:,it)=sim.Q(:,it-1)+sim.DQ(:,it);

    %Update active product lines
    sim.Iact(:,it)=(sim.Q(:,it)>eq.qmin);
    sim.AL(:,it)=sim.AL(:,it-1).*sim.Iact(:,it);

    % Growth
    sim.g(it)=mean(sim.Q(sim.Iact(:,it),it))/sim.Qbar(it-1)-1;

    % Remove growth from Q
    sim.Q(:,it)= sim.Q(:,it) /(1+sim.g(it));
    %sim.Q(:,it)= sim.Q(:,it) /(1+eq.g);

    %Mean productivity
    sim.Qbar(it)=mean(sim.Q(sim.Iact(:,it),it));

    %Update probability of internal innovation
    sim.x_a(:,it)=min(eq.x*dt*(sim.Q(:,it)/sim.Qbar(it)).^(par.alpha_x+par.gamma_x),1);
    sim.Px(:,it)= binornd(1,sim.x_a(:,it),[nprod,1]);

    %Update probability of entrant innovation
    sim.xe_a(:,it)=min(eq.xe*dt,1); %Adjusted so probability is bounded
    sim.Pe(:,it)= binornd(1,sim.xe_a(:,it),[nprod,1]);

    %p90 productivity
    sim.Qp90(it)=prctile(sim.Q(:,it),90);
    
    %Mean productivity Top 10
    sim.QbarT(it)=mean(sim.Q(sim.Q(:,it)>sim.Qp90(it),it));
    
    %Mean productivity Bottom 90
    sim.QbarB(it)=mean(sim.Q(sim.Q(:,it)<=sim.Qp90(it),it));


    % Stopping criteria
    if mod(it,dt_con)==0 && it>2*dt_con

           x=it-dt_con:it;
           y=sim.QbarT(x)./sim.Qbar(x);

       %Run regression
       sim.reg_con=fitlm(x,y);

       sim.res_conc=abs(table2array(sim.reg_con.Coefficients(2,1)))*dt_con;
       sim.res_conp=table2array(sim.reg_con.Coefficients(2,4));

    end

    
end


%Number of periods to converge
sim.Tsim=it;

%Test if any are NaN
if any(isnan(sim.Q(:,sim.Tsim)))

    disp('NaN in simulation');
    sim.Q(:,sim.Tsim);

end

    %----------------------------------------------------------------------
    %% Simulation variables
    %----------------------------------------------------------------------

    %Average growth
    sim.gbar=100*mean(sim.g(dt_res:end)/dt);

    %----------------------------------------------------------------------
    %% Save distribution
    %----------------------------------------------------------------------

    if save_Qdist==1
        Q0=sim.Q(:,sim.Tsim);
        save('sim_intq_B_Q0.mat','Q0');
    end

end