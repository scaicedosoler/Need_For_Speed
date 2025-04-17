function cf_ent_xq(par)

%Values of chi_b
nb=20;

chi_bv=linspace(0,10,nb);
chi_ev=par.chi_e*[0.6,1,2,3];

ne=length(chi_ev);

%Initialize vectors
xv=nan(nb,ne);
Qv=nan(nb,ne);
dqv=nan(nb,ne);
gv=nan(nb,ne);
xev=nan(nb,ne);
spibv=nan(nb,ne);

%Loop to solve the models

for ie=1:ne

    %Update entry parameter
    par.chi_e=chi_ev(ie);

    disp('---------------------------------------------')
    disp(['Running: chi_e=' num2str(par.chi_e) ])
    disp('---------------------------------------------')

    for ib=1:nb

        %Update benefit parameter
        par.chi_b=chi_bv(ib);

        %Run equilibrium
        eq=eq_sim_fun_intq_B(par);

        try
        %Fill vectors
        xv(ib,ie)=eq.x;
        Qv(ib,ie)=eq.Q;
        dqv(ib,ie)=eq.dq;
        gv(ib,ie)=eq.g;
        xev(ib,ie)=eq.xe;
        spibv(ib,ie)=100*eq.pib/(eq.pi+eq.pib+eq.piq);
        catch
            disp('Error')
        end

    end

end

%Compute percentage of entry
pe=round(100*(xev./(xev+xv)),0);


%Plots

%Speed
fig.cf_ent_xq_x=figure; hold all;
p1=plot(chi_bv, xv(:,1),'-'); p1.Color=par.opt.light_blue;
p2=plot(chi_bv, xv(:,2),'--'); p2.Color=par.opt.maroon;
p3=plot(chi_bv, xv(:,3),'-s'); p3.Color=par.opt.green;
p4=plot(chi_bv, xv(:,4),'-o'); p4.Color=par.opt.orange;
xlabel('Private Benefit $\chi_b$'); ylabel('Speed');
lgd=legend([p1,p2,p3,p4],{['Initial Entry ' num2str(pe(1,1)) '\%' ],['Initial Entry ' num2str(pe(1,2)) '\%' ],['Initial Entry ' num2str(pe(1,3)) '\%' ],['Initial Entry ' num2str(pe(1,4)) '\%']}); 
legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; 

%Quality
fig.cf_ent_xq_Q=figure; hold all;
p1=plot(chi_bv, Qv(:,1),'-'); p1.Color=par.opt.light_blue;
p2=plot(chi_bv, Qv(:,2),'--'); p2.Color=par.opt.maroon;
p3=plot(chi_bv, Qv(:,3),'-s'); p3.Color=par.opt.green;
p4=plot(chi_bv, Qv(:,4),'-o'); p4.Color=par.opt.orange;
xlabel('Private Benefit $\chi_b$'); ylabel('Quality');
lgd=legend([p1,p2,p3,p4],{['Initial Entry ' num2str(pe(1,1)) '\%' ],['Initial Entry ' num2str(pe(1,2)) '\%' ],['Initial Entry ' num2str(pe(1,3)) '\%' ],['Initial Entry ' num2str(pe(1,4)) '\%']});  
legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; 


%Growth
fig.cf_ent_xq_g=figure; hold all;
p1=plot(chi_bv, 100*gv(:,1),'-'); p1.Color=par.opt.light_blue;
p2=plot(chi_bv, 100*gv(:,2),'--'); p2.Color=par.opt.maroon;
p3=plot(chi_bv, 100*gv(:,3),'-s'); p3.Color=par.opt.green;
p4=plot(chi_bv, 100*gv(:,4),'-o'); p4.Color=par.opt.orange;
xlabel('Private Benefit $\chi_b$'); ylabel('Growth \%');
lgd=legend([p1,p2,p3,p4],{['Initial Entry ' num2str(pe(1,1)) '\%' ],['Initial Entry ' num2str(pe(1,2)) '\%' ],['Initial Entry ' num2str(pe(1,3)) '\%' ],['Initial Entry ' num2str(pe(1,4)) '\%']}); 
legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; 

%-----------------
% Relative to share of profits
%-----------------

%Speed
fig.cf_ent_xq_x_spib=figure; hold all;
p1=plot(spibv(:,1), xv(:,1),'-'); p1.Color=par.opt.light_blue;
p2=plot(spibv(:,2), xv(:,2),'--'); p2.Color=par.opt.maroon;
p3=plot(spibv(:,3), xv(:,3),'-s'); p3.Color=par.opt.green;
p4=plot(spibv(:,4), xv(:,4),'-o'); p4.Color=par.opt.orange;
xlabel('Share of Profits of Private Benefit (\%)'); ylabel('Speed');
lgd=legend([p1,p2,p3,p4],{['Initial Entry ' num2str(pe(1,1)) '\%' ],['Initial Entry ' num2str(pe(1,2)) '\%' ],['Initial Entry ' num2str(pe(1,3)) '\%' ],['Initial Entry ' num2str(pe(1,4)) '\%']}); 
legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; 

%Quality
fig.cf_ent_xq_Q_spib=figure; hold all;
p1=plot(spibv(:,1), Qv(:,1),'-'); p1.Color=par.opt.light_blue;
p2=plot(spibv(:,2), Qv(:,2),'--'); p2.Color=par.opt.maroon;
p3=plot(spibv(:,3), Qv(:,3),'-s'); p3.Color=par.opt.green;
p4=plot(spibv(:,4), Qv(:,4),'-o'); p4.Color=par.opt.orange;
xlabel('Share of Profits of Private Benefit (\%)'); ylabel('Quality');
lgd=legend([p1,p2,p3,p4],{['Initial Entry ' num2str(pe(1,1)) '\%' ],['Initial Entry ' num2str(pe(1,2)) '\%' ],['Initial Entry ' num2str(pe(1,3)) '\%' ],['Initial Entry ' num2str(pe(1,4)) '\%']});  
legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; 


%Growth
fig.cf_ent_xq_g_spib=figure; hold all;
p1=plot(spibv(:,1), 100*gv(:,1),'-'); p1.Color=par.opt.light_blue;
p2=plot(spibv(:,2), 100*gv(:,2),'--'); p2.Color=par.opt.maroon;
p3=plot(spibv(:,3), 100*gv(:,3),'-s'); p3.Color=par.opt.green;
p4=plot(spibv(:,4), 100*gv(:,4),'-o'); p4.Color=par.opt.orange;
xlabel('Share of Profits of Private Benefit (\%)'); ylabel('Growth \%');
lgd=legend([p1,p2,p3,p4],{['Initial Entry ' num2str(pe(1,1)) '\%' ],['Initial Entry ' num2str(pe(1,2)) '\%' ],['Initial Entry ' num2str(pe(1,3)) '\%' ],['Initial Entry ' num2str(pe(1,4)) '\%']}); 
legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; 

%% Save figures
%par.opt.save_fig=1;
if ~isfield(par,'save_str')
    par.save_str='';
end

if par.opt.save_fig==1
    exportgraphics(fig.cf_ent_xq_x,[par.opt.dir_fig 'cf_ent_xq_x' par.save_str '.pdf'])
    exportgraphics(fig.cf_ent_xq_Q,[par.opt.dir_fig 'cf_ent_xq_Q' par.save_str '.pdf'])
    exportgraphics(fig.cf_ent_xq_g,[par.opt.dir_fig 'cf_ent_xq_g' par.save_str '.pdf'])
    
    exportgraphics(fig.cf_ent_xq_x_spib,[par.opt.dir_fig 'cf_ent_xq_x_spib' par.save_str '.pdf'])
    exportgraphics(fig.cf_ent_xq_Q_spib,[par.opt.dir_fig 'cf_ent_xq_Q_spib' par.save_str '.pdf'])
    exportgraphics(fig.cf_ent_xq_g_spib,[par.opt.dir_fig 'cf_ent_xq_g_spib' par.save_str '.pdf'])
end

