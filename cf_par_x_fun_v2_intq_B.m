%Computes the counterfactual of parameters depending on x
% Version 2: How much of the change can be attributed to the spillover
% For the intq_B version of the model
function cf_par_x_fun_v2_intq_B(par)

%Inputs
% par: baseline parameters
% eq1: The comparison equilibrium structure

%Change in Speed (from data)
dx=1.5553;

%Run for baseline parameters
eq=eq_sim_fun_intq_B(par);

%% Change of chi_b

%Base parameters
parb=par;

%Counterfactual with chi_bc
%par.chi_b=0;

%Parameters for counterfactuals
ng=9;

zeta_est=0.12;
zeta_min=0; zeta_max=0.25; 

zetag=sort([linspace(zeta_min,zeta_max,ng),zeta_est]);

ind_est=find(zetag==zeta_est);

ng=length(zetag);


%% Variables

%Initialize vectors
gc=NaN(4,ng);
slec=NaN(4,ng);
sec=NaN(4,ng);
q_ent_incc=NaN(4,ng);
sl_topc=NaN(4,ng);
s_topc=NaN(4,ng);
q_top_botc=NaN(4,ng);
thetagc=NaN(4,ng);


%% Compute moments for each parameter combination

for ig=1:ng

    % Quality: lambda
    par=parb;
    %par.chi_b=0;
    par.lambda=parb.lambda*dx^(-zetag(ig));    
    eqc=eq_sim_fun_intq_B(par);
    gc(1,ig)=eqc.g;
    slec(1,ig)=eqc.sle;
    sec(1,ig)=eqc.se;
    q_ent_incc(1,ig)=eqc.q_ent_inc;
    sl_topc(1,ig)=eqc.sl_top;
    s_topc(1,ig)=eqc.s_top;
    q_top_botc(1,ig)=eqc.q_top_bot;
    thetagc(1,ig)=par.lambda;


    % Entrant productivity: chi_e
    par=parb;
    %par.chi_b=0;
    par.chi_e=parb.chi_e*dx^(-zetag(ig));    
    eqc=eq_sim_fun_intq_B(par);
    gc(2,ig)=eqc.g;
    slec(2,ig)=eqc.sle;
    sec(2,ig)=eqc.se;
    q_ent_incc(2,ig)=eqc.q_ent_inc;
    sl_topc(2,ig)=eqc.sl_top;
    s_topc(2,ig)=eqc.s_top;
    q_top_botc(2,ig)=eqc.q_top_bot;
    thetagc(2,ig)=par.chi_e;



    % Entrant quality: lambda_e
    par=parb;
    %par.chi_b=0;
    par.lambda_e=parb.lambda_e*dx^(-zetag(ig));   
    eqc=eq_sim_fun_intq_B(par);
    gc(3,ig)=eqc.g;
    slec(3,ig)=eqc.sle;
    sec(3,ig)=eqc.se;
    q_ent_incc(3,ig)=eqc.q_ent_inc;
    sl_topc(3,ig)=eqc.sl_top;
    s_topc(3,ig)=eqc.s_top;
    q_top_botc(3,ig)=eqc.q_top_bot;
    thetagc(3,ig)=par.lambda_e;

    % Knowledge diffusion: nu
    par=parb;
    %par.chi_b=0;
    par.nu=parb.nu*dx^(-zetag(ig));    
    eqc=eq_sim_fun_intq_B(par);
    gc(4,ig)=eqc.g;
    slec(4,ig)=eqc.sle;
    sec(4,ig)=eqc.se;
    q_ent_incc(4,ig)=eqc.q_ent_inc;
    sl_topc(4,ig)=eqc.sl_top;
    s_topc(4,ig)=eqc.s_top;
    q_top_botc(4,ig)=eqc.q_top_bot;
    thetagc(4,ig)=par.nu;

end

%% Plot
ng=length(zetag);

%Growth
fig_g=figure; hold all;
p0=plot(zetag,100*eq.g*ones(1,ng),':k');
p1=plot(zetag,100*gc(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),100*gc(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,100*gc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,100*gc(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Growth (\%)');
%xlim([0,0.5])

%Entrant Labor
fig_sle=figure; hold all;
p0=plot(zetag,100*eq.sle*ones(1,ng),':k');
p1=plot(zetag,100*slec(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),100*slec(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,100*slec(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,100*slec(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Labor for Entrants (\%)');
legend("off")
%xlim([0,0.5])

%Entrant Innovation 
fig_se=figure; hold all;
p0=plot(zetag,100*eq.se*ones(1,ng),':k');
p1=plot(zetag,100*sec(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),100*sec(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,100*sec(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,100*sec(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Entrant Innovation (\%)');
legend("off")
%xlim([0,0.5])

%Entrant Quality versus Incumbent Quality
fig_q_ent_inc=figure; hold all;
p0=plot(zetag,eq.q_ent_inc*ones(1,ng),':k');
p1=plot(zetag,q_ent_incc(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),q_ent_incc(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,q_ent_incc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,q_ent_incc(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Entrant/Incumbent Quality ');
legend("off")
%xlim([0,0.5])

%Labor Concentration in Top 10
fig_sl_top=figure; hold all;
p0=plot(zetag,100*eq.sl_top*ones(1,ng),':k');
p1=plot(zetag,100*sl_topc(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),100*sl_topc(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,100*sl_topc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,100*sl_topc(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Share of Labor in Top 10 (\%)');
legend("off")
%xlim([0,0.5])

%Innovation Concentration in Top 10
fig_s_top=figure; hold all;
p0=plot(zetag,100*eq.s_top*ones(1,ng),':k');
p1=plot(zetag,100*s_topc(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),100*s_topc(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,100*s_topc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,100*s_topc(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Share of Innovations in Top 10 (\%)');
legend("off")
%xlim([0,0.5])

%Quality Top 10/ Bottom 90
fig_q_top_bot=figure; hold all;
p0=plot(zetag,eq.q_top_bot*ones(1,ng),':k');
p1=plot(zetag,q_top_botc(1,:),'-'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
s1=scatter(zetag(ind_est),q_top_botc(1,ind_est)); s1.Marker="square";s1.MarkerFaceColor='w'; s1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(zetag,q_top_botc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(zetag,q_top_botc(3,:),'-.'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Elasticity $\zeta$'); ylabel('Top 10/ Bottom 90 Quality');
legend("off")
%xlim([0,0.5])


%% Save plots

if par.opt.save_fig==1

    %For Matlab 2020 onward

    exportgraphics(fig_g,[par.opt.dir_fig  'cf_par_x_v2_g_' par.opt.save_str '.pdf'])
    exportgraphics(fig_sle,[par.opt.dir_fig  'cf_par_x_v2_sle_' par.opt.save_str '.pdf'])
    exportgraphics(fig_se,[par.opt.dir_fig  'cf_par_x_v2_se_' par.opt.save_str '.pdf'])
    exportgraphics(fig_q_ent_inc,[par.opt.dir_fig  'cf_par_x_v2_q_ent_inc_' par.opt.save_str '.pdf'])
    exportgraphics(fig_sl_top,[par.opt.dir_fig  'cf_par_x_v2_sl_top_' par.opt.save_str '.pdf'])
    exportgraphics(fig_s_top,[par.opt.dir_fig  'cf_par_x_v2_s_top_' par.opt.save_str '.pdf'])
    exportgraphics(fig_q_top_bot,[par.opt.dir_fig  'cf_par_x_v2_q_top_bot_' par.opt.save_str '.pdf'])

    %close all;

end

end