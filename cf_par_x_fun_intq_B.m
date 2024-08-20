%Computes the counterfactual of parameters depending on x
% For the intq_B version of the model
function cf_par_x_fun_intq_B(par)

%Run for baseline parameters
eq=eq_sim_fun_intq_B(par);

%% Change of chi_b

%Base parameters
parb=par;

%Counterfactual with chi_bc
par.chi_b=0;

%Parameters for counterfactuals
ng=10;
ratio_min=1; ratio_max=2; 
ratiog=linspace(ratio_min,ratio_max,ng);

%Compute the equivalent zeta
zetag=log(ratiog)/eq.x;


%% Variables

%Initialize vectors
gc=NaN(4,ng);
slec=NaN(4,ng);
sec=NaN(4,ng);
q_ent_incc=NaN(4,ng);
sl_topc=NaN(4,ng);
s_topc=NaN(4,ng);
q_top_botc=NaN(4,ng);




%% Compute moments for each parameter combination

for ig=1:ng

    % Quality: lambda
    par=parb;
    par.chi_b=0;
    par.lambda=parb.lambda*ratiog(ig);    
    eqc=eq_sim_fun_intq_B(par);
    gc(1,ig)=eqc.g;
    slec(1,ig)=eqc.sle;
    sec(1,ig)=eqc.se;
    q_ent_incc(1,ig)=eqc.q_ent_inc;
    sl_topc(1,ig)=eqc.sl_top;
    s_topc(1,ig)=eqc.s_top;
    q_top_botc(1,ig)=eqc.q_top_bot;


    % Entrant productivity: chi_e
    par=parb;
    par.chi_b=0;
    par.chi_e=parb.chi_e*ratiog(ig);    
    eqc=eq_sim_fun_intq_B(par);
    gc(2,ig)=eqc.g;
    slec(2,ig)=eqc.sle;
    sec(2,ig)=eqc.se;
    q_ent_incc(2,ig)=eqc.q_ent_inc;
    sl_topc(2,ig)=eqc.sl_top;
    s_topc(2,ig)=eqc.s_top;
    q_top_botc(2,ig)=eqc.q_top_bot;


    % Entrant quality: lambda_e
    par=parb;
    par.chi_b=0;
    par.lambda_e=parb.lambda_e*ratiog(ig);    
    eqc=eq_sim_fun_intq_B(par);
    gc(3,ig)=eqc.g;
    slec(3,ig)=eqc.sle;
    sec(3,ig)=eqc.se;
    q_ent_incc(3,ig)=eqc.q_ent_inc;
    sl_topc(3,ig)=eqc.sl_top;
    s_topc(3,ig)=eqc.s_top;
    q_top_botc(3,ig)=eqc.q_top_bot;

    % Knowledge diffusion: nu
    par=parb;
    par.chi_b=0;
    par.nu=parb.nu*ratiog(ig);    
    eqc=eq_sim_fun_intq_B(par);
    gc(4,ig)=eqc.g;
    slec(4,ig)=eqc.sle;
    sec(4,ig)=eqc.se;
    q_ent_incc(4,ig)=eqc.q_ent_inc;
    sl_topc(4,ig)=eqc.sl_top;
    s_topc(4,ig)=eqc.s_top;
    q_top_botc(4,ig)=eqc.q_top_bot;

end

%% Plot
dprcg=100*(ratiog-1);


%Growth
fig_g=figure; hold all;
p0=plot(dprcg,100*eq.g*ones(1,ng),':k');
p1=plot(dprcg,100*gc(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,100*gc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,100*gc(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Growth (\%)');

%Entrant Labor
fig_sle=figure; hold all;
p0=plot(dprcg,100*eq.sle*ones(1,ng),':k');
p1=plot(dprcg,100*slec(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,100*slec(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,100*slec(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Labor for Entrants (\%)');
legend("off")

%Entrant Innovation 
fig_se=figure; hold all;
p0=plot(dprcg,100*eq.se*ones(1,ng),':k');
p1=plot(dprcg,100*sec(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,100*sec(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,100*sec(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Entrant Innovation (\%)');
legend("off")

%Entrant Quality versus Incumbent Quality
fig_q_ent_inc=figure; hold all;
p0=plot(dprcg,eq.q_ent_inc*ones(1,ng),':k');
p1=plot(dprcg,q_ent_incc(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,q_ent_incc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,q_ent_incc(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Entrant/Incumbent Quality ');
legend("off")

%Labor Concentration in Top 10
fig_sl_top=figure; hold all;
p0=plot(dprcg,100*eq.sl_top*ones(1,ng),':k');
p1=plot(dprcg,100*sl_topc(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,100*sl_topc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,100*sl_topc(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Share of Labor in Top 10 (\%)');
legend("off")

%Innovation Concentration in Top 10
fig_s_top=figure; hold all;
p0=plot(dprcg,100*eq.s_top*ones(1,ng),':k');
p1=plot(dprcg,100*s_topc(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,100*s_topc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,100*s_topc(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Share of Innovations in Top 10 (\%)');
legend("off")

%Quality Top 10/ Bottom 90
fig_q_top_bot=figure; hold all;
p0=plot(dprcg,eq.q_top_bot*ones(1,ng),':k');
p1=plot(dprcg,q_top_botc(1,:),'-O'); p1.Color=par.opt.light_blue; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par.opt.light_blue;
p2=plot(dprcg,q_top_botc(2,:),'--'); p2.Color=par.opt.maroon;
p3=plot(dprcg,q_top_botc(3,:),'-D'); p3.Color=par.opt.green; p3.MarkerFaceColor='w'; p3.MarkerEdgeColor=par.opt.green;
lgd=legend([p0,p1,p2, p3],{['Calibration ' num2str(par.iyear) '-' num2str(par.iyear+5)],'Quality Level $\lambda$','Entrant Productivity $\chi_e$','Entrant Quality $\lambda_e$'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Percentage Change in Parameter (\%)'); ylabel('Top 10/ Bottom 90 Quality');
legend("off")


%% Save plots

if par.opt.save_fig==1

    %For Matlab 2020 onward

    exportgraphics(fig_g,[par.opt.dir_fig  'cf_par_x_g_' par.opt.save_str '.pdf'])
    exportgraphics(fig_sle,[par.opt.dir_fig  'cf_par_x_sle_' par.opt.save_str '.pdf'])
    exportgraphics(fig_se,[par.opt.dir_fig  'cf_par_x_se_' par.opt.save_str '.pdf'])
    exportgraphics(fig_q_ent_inc,[par.opt.dir_fig  'cf_par_x_q_ent_inc_' par.opt.save_str '.pdf'])
    exportgraphics(fig_sl_top,[par.opt.dir_fig  'cf_par_x_sl_top_' par.opt.save_str '.pdf'])
    exportgraphics(fig_s_top,[par.opt.dir_fig  'cf_par_x_s_top_' par.opt.save_str '.pdf'])
    exportgraphics(fig_q_top_bot,[par.opt.dir_fig  'cf_par_x_q_top_bot_' par.opt.save_str '.pdf'])

    %close all;

end


end