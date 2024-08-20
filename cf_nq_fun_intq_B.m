%Computes the counterfactual of no endogenous quality
% For the intq_B version of the model
function cf_nq_fun_intq_B(smm)

%Unpack structure and parameters
par=smm.par;

%No simulations
par.opt.run_sim=0;

%Grid for chi_b
chi_bg=0:0.2:max(10,par.chi_b);
chi_bg=sort(unique([chi_bg, par.chi_b]));

ng=length(chi_bg);

%% Moments

%Initialize vectors
g1=NaN(1,ng);gc=NaN(1,ng);
spib1=NaN(1,ng);spibc=NaN(1,ng);
slx1=NaN(1,ng);slxc=NaN(1,ng);

for ig=1:ng

    try
        %Baseline
        par.chi_b=chi_bg(ig);
        eq1=eq_sim_fun_intq_B(par);

        %No Endogenous quality
        parc=par;
        parc.lambda=par.lambda_nq;

        eqc=eq_sim_fun_intq_B_nq(parc);

        %Vairables
        g1(ig)=eq1.g;
        spib1(ig)=eq1.pib/(eq1.pi+eq1.pib+eq1.piq);
        slx1(ig)=eq1.lx/par.L_I;

        gc(ig)=eqc.g;
        spibc(ig)=eqc.pib/(eqc.pi+eqc.pib+eqc.piq);
        slxc(ig)=eqc.lx/parc.L_I;

    catch
        disp('No convergence')
    end

end



%% Plots

posp=isfinite(g1);

%Growth
hfig.g=figure; hold all;
p1=plot(chi_bg(posp),100*g1(posp)); p1.Color=par.opt.light_blue;
pc=plot(chi_bg(posp),100*gc(posp),'--'); pc.Color=par.opt.maroon;
lgd=legend([p1,pc],{'Endogenous Quality','Fixed Quality'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Private Benefit $\chi_b$'); ylabel('Growth (\%)');

%Share of profits
hfig.spib=figure; hold all;
p1=plot(chi_bg(posp),100*spib1(posp)); p1.Color=par.opt.light_blue;
pc=plot(chi_bg(posp),100*spibc(posp),'--'); pc.Color=par.opt.maroon;
lgd=legend([p1,pc],{'Endogenous Quality','Fixed Quality'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Private Benefit $\chi_b$'); ylabel('Share of Profits of Private Benefit (\%)');

%Share of labor
hfig.slx=figure; hold all;
p1=plot(chi_bg(posp),100*slx1(posp)); p1.Color=par.opt.light_blue;
pc=plot(chi_bg(posp),100*slxc(posp),'--'); pc.Color=par.opt.maroon;
lgd=legend([p1,pc],{'Endogenous Quality','Fixed Quality'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Private Benefit $\chi_b$'); ylabel('Share of Labor for Speed (\%)');

%Growth vs share of profits
hfig.spib_g=figure; hold all;
p1=plot(100*spib1(posp),100*g1(posp)); p1.Color=par.opt.light_blue;
pc=plot(100*spibc(posp),100*gc(posp),'--'); pc.Color=par.opt.maroon;
lgd=legend([p1,pc],{'Endogenous Quality','Fixed Quality'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
xlim([0,min(100*max(spib1),100*max(spibc))])
xlabel('Share of Profits of Private Benefit (\%)'); ylabel('Growth (\%)');



%% Save plots

if par.opt.save_fig==1

    %For Matlab 2020 onward

    exportgraphics(hfig.g,[par.opt.dir_fig  'cf_nq_g_intq_B' par.save_str '.pdf'])
    exportgraphics(hfig.spib,[par.opt.dir_fig 'cf_nq_spib_intq_B' par.save_str '.pdf'])
    exportgraphics(hfig.slx,[par.opt.dir_fig  'cf_nq_slx_intq_B' par.save_str '.pdf'])
    exportgraphics(hfig.spib_g,[par.opt.dir_fig  'cf_nq_spib_g_intq_B' par.save_str '.pdf'])

    %close all;

end


end