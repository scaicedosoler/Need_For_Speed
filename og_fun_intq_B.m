%Computes the optimal growth and the equilibrium of a given parameter
% For the intq_B version of the model
function og=og_fun_intq_B(par,og)

%Options
par.opt.run_sim=0; %Do not run simulations

%Unpack og parameter
parg=og.([og.par 'g']);

if  isfield(og,'greek') && og.greek==0
    par_lbl=['$' og.par '$'];
    if og.dsub==1
        par_lbl=[ '$' insertAfter(og.par,'_','{') '}$'];
    end
else
    par_lbl=['$\' og.par '$'];
    if og.dsub==1
        par_lbl=[ '$\' insertAfter(og.par,'_','{') '}$'];
    end
end


%Save initial value and position
par0=par.(og.par);
pos0=find(parg==par0);

%Initialize variables
ng=length(parg);

%Variables
%Equilibrium
gg=NaN(1,ng); slxg=NaN(1,ng); slqg=NaN(1,ng);  sleg=NaN(1,ng); xg=NaN(1,ng); Qg=NaN(1,ng); xeg=NaN(1,ng);

%Social planner
gog=NaN(1,ng); slxog=NaN(1,ng); slqog=NaN(1,ng);  sleog=NaN(1,ng); xog=NaN(1,ng); Qog=NaN(1,ng); xeog=NaN(1,ng);

%Comparison
r_gg=NaN(1,ng);


%Diagnostic parameters
noconv_par_eq=NaN(1,ng);
noconv_par_sp=NaN(1,ng);

for ig=1:ng
    par.(og.par)=parg(ig);

    %------------------
    %Compute equilibrium
    %------------------

    try
        eq=eq_sim_fun_intq_B(par);

        gg(ig)=eq.g;
        slxg(ig)=eq.lx/par.L_I;
        slqg(ig)=eq.lq/par.L_I;
        sleg(ig)=eq.le/par.L_I;
        xg(ig)=eq.x;
        Qg(ig)=eq.dq;
        xeg(ig)=eq.xe;

    catch
        %Not converging parameters
        disp('Parameter for equilibiurm not converging')
        fprintf('-------------------------\n')
        noconv_par_eq(ig)=parg(ig);
    end

    %----------------------
    %Compute social planner
    %----------------------

    try

        %Compute equilibrium
        sp=sp_intq_B(par);

        gog(ig)=sp.g;
        slxog(ig)=sp.lx/par.L_I;
        slqog(ig)=sp.lq/par.L_I;
        sleog(ig)=sp.le/par.L_I;
        xog(ig)=sp.x;
        Qog(ig)=sp.Q;
        xeog(ig)=sp.xe;

    catch
        %Not converging parameters
        disp('Parameter for social planner not converging')
        fprintf('-------------------------\n')
        noconv_par_sp(ig)=parg(ig);
    end

    %Comparison
    r_gg(ig)=gg(ig)/gog(ig);

end


%Save all in cell
og.var={gg,slxg,slqg,sleg,xg,Qg,xeg,...
        gog,slxog,slqog,sleog,xog,Qog,xeog...
        r_gg};


%Save just in case
save([par.opt.dir_mat 'og_' og.par '.mat'])


%% Plots
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
ind=isfinite(gg);
indo=isfinite(gog);

%Moments
hfig.g=figure; hold all;
p=plot(parg(ind),gg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),gg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),gog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),gog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Growth $g$');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';


hfig.slx=figure; hold all;
p=plot(parg(ind),slxg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),slxg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),slxog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),slxog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Share of Labor for Speed');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';

hfig.slq=figure; hold all;
p=plot(parg(ind),slqg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),slqg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),slqog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),slqog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Share of Labor for Quality');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';

hfig.sle=figure; hold all;
p=plot(parg(ind),sleg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),sleg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),sleog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),sleog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Share of Labor for Entrants');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';

hfig.x=figure; hold all;
p=plot(parg(ind),xg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),xg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),xog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),xog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Arrival Rate');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';

hfig.Q=figure; hold all;
p=plot(parg(ind),Qg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),Qg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),Qog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),Qog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Quality');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';

hfig.xe=figure; hold all;
p=plot(parg(ind),xeg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),xeg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
po=plot(parg(indo),xeog(indo),'--'); po.Color=par.opt.maroon;
so1=scatter(parg(pos0),xeog(pos0),'D'); so1.MarkerFaceColor=par.opt.maroon; so1.MarkerEdgeColor='w';
%ylim([0.99*min(gg),1.01*max(gg)])
xlabel(par_lbl); ylabel('Entrants Arrival Rate');
lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';

hfig.r_g=figure; hold all;
p=plot(parg(ind),r_gg(ind)); p.Color=par.opt.light_blue;
s1=scatter(parg(pos0),r_gg(pos0),'o'); s1.MarkerFaceColor=par.opt.light_blue; s1.MarkerEdgeColor='w';
xlabel(par_lbl); ylabel('Growth Rates Ratio: $g^*/g^s$');
%lgd=legend([p,po],{'Competitive Equilibrium','Optimal Growth'});lgd.EdgeColor="none"; lgd.Location='best';


%% Growth relative to chi_b

%Grid for chi_b
nchib=20; chib_min=0; chib_max=40;
chibg=linspace(chib_min,chib_max,nchib);

%Growth rates ratio
hfig.g_chib=figure; hold all; xlabel('Private Benefit $\chi_b$'); ylabel('Growth Rates Ratio: $g^*/g^s$');

%Use the succesful parameters
J=find(~isnan(gg));

%Limit to 5 lines
if length(J)>4
    J=unique([min(J), ceil(max(J)/4), ceil(max(J)/2), ceil(3*max(J)/4), max(J)]);
end

%lstyle={'-.','--',':','-.'} ;
lstyle={'-','--','-o','-D','-s'} ;
lcolor={par.opt.light_blue,par.opt.maroon, par.opt.green, 'black', par.opt.orange};

pv1=[]; 



for jg=1:length(J)

    %Upadate parameter
     par.(og.par)=parg(jg);

     %Initialize vectors
        ggc=NaN(1,nchib); gogc=NaN(1,nchib); r_ggc=NaN(1,nchib);

    for ichib=1:nchib
        %Update chi_b
        par.chi_b=chibg(ichib);

        try
            eqc=eq_sim_fun_intq_B(par);
            ggc(ichib)=eqc.g;

            spc=sp_intq_B(par);
            gogc(ichib)=spc.g;

        catch
            %Not converging parameters
            disp('Parameter for equilibiurm or social planner not converging')
            fprintf('-------------------------\n')
        end

        %Ratio
        r_ggc(ichib)=ggc(ichib)/gogc(ichib);

    end

    %Plot moment functions
    figure(hfig.g_chib);
    p1(jg)=plot(chibg,r_ggc,lstyle{jg}); p1(jg).Color=lcolor{jg}; p1(jg).MarkerFaceColor='w';

    %For legend
    pv1=[pv1, p1(jg)];  lgd_c{jg}=[par_lbl '=' num2str(round(parg(J(jg)),2))];
end

%Legend
figure(hfig.g_chib);
lgd=legend(pv1, lgd_c);
lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
ylim([0.8,1])


%% Save plots
%par.opt.save_fig=1;

if ~isfield(par.opt,'save_str2')
    par.opt.save_str2='';
end

if par.opt.save_fig==1

    %For Matlab 2020 onward
    exportgraphics(hfig.g,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_g_intq_B.pdf'])
    exportgraphics(hfig.slx,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_slx_intq_B.pdf'])
    exportgraphics(hfig.slq,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_slq_intq_B.pdf'])
    exportgraphics(hfig.sle,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_sle_intq_B.pdf'])
    exportgraphics(hfig.x,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_x_intq_B.pdf'])
    exportgraphics(hfig.Q,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_Q_intq_B.pdf'])
    exportgraphics(hfig.xe,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_xe_intq_B.pdf'])

    exportgraphics(hfig.g_chib,[par.opt.dir_fig 'og_' par.opt.save_str2 og.par '_g_chib_intq_B.pdf'])

    %close all;

end


end