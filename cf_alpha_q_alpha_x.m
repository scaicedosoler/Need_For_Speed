%Counterfactual: Changing alpha_q and alpha_x
function cf_alpha_q_alpha_x(par,eq,par0,eq0)

if ~isfield(par,'byear')
    par.byear=2010;
    par.cyear=1980;
end

%Choose base parameters
parb=par;

%Number of parameters
npar=3;

%Vector of parameters
% alpha_qv=sort([linspace(0.01,0.6,npar),par.alpha_q]);
% alpha_xv=sort([linspace(0.1,0.6,npar),par.alpha_x]);

alpha_qv=[par.alpha_q,0.01,0.31,0.6,];
alpha_xv=[par.alpha_x,0.01,0.08,0.6,];

npar=length(alpha_qv);

%Number of points for plots
ng=999;

%Incumbent labor
linc0=eq0.lx+eq0.lq;
linc=eq.lx+eq.lq;

lxg=sort([linspace(0.01*linc,0.99*linc,ng),eq.lx]);
lx0g=sort([linspace(0.01*linc0,0.99*linc0,ng),eq0.lx]);

ind=find(lxg==eq.lx);
ind0=find(lx0g==eq0.lx);

plxg=lxg./linc;
plx0g=lx0g./linc0;

lqg=(1-plxg)*linc;
lq0g=(1-plx0g)*linc0;

% %Arrival rate
% xg=par.chi*lxg.^par.alpha_x;
% x0g=par0.chi*lx0g.^par0.alpha_x;
%
% %Quality
% dqg=par.lambda*lqg.^par.alpha_q;
% dq0g=par0.lambda*lq0g.^par0.alpha_q;
%
% %Entry
% xe=par.chi_e*eq.le.^par.alpha_e;
% dqe=par.lambda_e;
%
% xe0=par0.chi_e*eq0.le.^par0.alpha_e;
% dqe0=par0.lambda_e;
%
% %Growth
% gg=xe.*dqe+xg.*dqg;
% g0g=xe0.*dqe0+x0g.*dq0g;

%% -------1. Quality elasticity of labor: alpha_q----------

%Initialize plots

%Speed
fig.xclxq_tot_alpha_q=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Arrival Rate')


%Quality
fig.dqclxq_tot_alpha_q=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Quality')


%Growth
fig.gclxq_tot_alpha_q=figure; hold all
%sb=scatter(eq.lx/par.L_I,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Growth \%')

%Normalized plots

%Speed
fig.xaclxq_tot_alpha_q=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Arrival Rate')


%Quality
fig.dqaclxq_tot_alpha_q=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Quality')


%Growth
fig.gaclxq_tot_alpha_q=figure; hold all
%sb=scatter(eq.lx/par.L_I,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Growth \%')

% Cell with colors
pcolorc={'k',par.opt.light_blue,par.opt.maroon,par.opt.green};
stylec={':','-','--','-.'};

lgds=[]; lgdq=[]; lgdg=[];
lgdsa=[]; lgdqa=[]; lgdga=[];
lgdstr=cell(1, npar);

%Calibration value
xcal=eq.x;
dqcal=eq.dq;
gcal=eq.g;


for ipar=1:npar

    %Update parameter
    par.alpha_q=alpha_qv(ipar);

    %Changing the speed and quality of labor (Relative to total labor)

    %Arrival rate
    xg=par.chi*lxg.^par.alpha_x;
    x0g=par0.chi*lx0g.^par0.alpha_x;

    %Quality
    dqg=par.lambda*lqg.^par.alpha_q;
    dq0g=par0.lambda*lq0g.^par0.alpha_q;

    %Entry
    xe=par.chi_e*eq.le.^par.alpha_e;
    dqe=par.lambda_e;

    xe0=par0.chi_e*eq0.le.^par0.alpha_e;
    dqe0=par0.lambda_e;

    %Growth
    gg=xe.*dqe+xg.*dqg;
    g0g=xe0.*dqe0+x0g.*dq0g;


    %Proportion relative to total

    %-----------------------------
    %Reallocation of Labor
    %-----------------------------

    %Speed
    figure(fig.xclxq_tot_alpha_q);
    ps(ipar)=plot(lxg/par.L_I,xg,stylec{ipar}); ps(ipar).Color=pcolorc{ipar};

    %Quality
    figure(fig.dqclxq_tot_alpha_q)
    pq(ipar)=plot(lxg/par.L_I,dqg,stylec{ipar}); pq(ipar).Color=pcolorc{ipar};


    %Growth
    figure(fig.gclxq_tot_alpha_q)
    pg(ipar)=plot(lxg/par.L_I,100*gg,stylec{ipar}); pg(ipar).Color=pcolorc{ipar};

    %Legends
    lgdstr{ipar}=['$\alpha_q$=' num2str(par.alpha_q,2)  ];
    lgds=[lgds,ps(ipar)]; lgdq=[lgdq,pq(ipar)];  lgdg=[lgdg,pg(ipar)];

    %-----------------------------
    %Normalizing to calibration
    %-----------------------------

    %Additive normalization
    xag=xg+(xcal-xg(ind));
    dqag=dqg+(dqcal-dqg(ind));
    gag=gg+(gcal-gg(ind));


    %Find max
    [~,indmax]=max(gag);
    plot_max=1;

    %Speed
    figure(fig.xaclxq_tot_alpha_q);
    psa(ipar)=plot(lxg/par.L_I,xag,stylec{ipar}); psa(ipar).Color=pcolorc{ipar};

    %Quality
    figure(fig.dqaclxq_tot_alpha_q)
    pqa(ipar)=plot(lxg/par.L_I,dqag,stylec{ipar}); pqa(ipar).Color=pcolorc{ipar};


    %Growth
    figure(fig.gaclxq_tot_alpha_q)
    pga(ipar)=plot(lxg/par.L_I,100*gag,stylec{ipar}); pga(ipar).Color=pcolorc{ipar};

    if plot_max==1
        plxg=lxg/par.L_I;
        smax=scatter(plxg(indmax),100*gag(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=pcolorc{ipar}; smax.MarkerEdgeColor=pcolorc{ipar};
    end


    %Legends
    if ipar==1
        lgdastr{ipar}=['Calibration $\alpha_q$=' num2str(par.alpha_q,2)  ];
        lgdsa=[lgdsa,psa(ipar)]; lgdqa=[lgdqa,pqa(ipar)];  lgdga=[lgdga,pga(ipar)];
    else
        lgdastr{ipar}=['$\alpha_q$=' num2str(par.alpha_q,2)  ];
        lgdsa=[lgdsa,psa(ipar)]; lgdqa=[lgdqa,pqa(ipar)];  lgdga=[lgdga,pga(ipar)];
    end


end

%Calibration lines and legend

figure(fig.xclxq_tot_alpha_q);
yLimits = get(gca,'YLim');
cx=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgds,lgdstr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="northwest";

figure(fig.dqclxq_tot_alpha_q)
yLimits = get(gca,'YLim');
cq=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdq,lgdstr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

figure(fig.gclxq_tot_alpha_q)
yLimits = get(gca,'YLim');
cg=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdg,lgdstr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

%Normalized
figure(fig.xaclxq_tot_alpha_q);
yLimits = get(gca,'YLim');
cx=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdsa,lgdastr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

figure(fig.dqaclxq_tot_alpha_q)
yLimits = get(gca,'YLim');
cq=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdqa,lgdastr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

figure(fig.gaclxq_tot_alpha_q)
yLimits = get(gca,'YLim');
cg=plot([eq.lx/par.L_I,eq.lx/par.L_I],[0.5,2],':k');
lgd=legend(lgdga,lgdastr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";
ylim([0.5,2])

%% -------2. Speed elasticity of labor: alpha_x------------

%Change parameters to original value
par=parb;

%Initialize plots

%Speed
fig.xclxq_tot_alpha_x=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Arrival Rate')


%Quality
fig.dqclxq_tot_alpha_x=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Quality')


%Growth
fig.gclxq_tot_alpha_x=figure; hold all
%sb=scatter(eq.lx/par.L_I,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Growth \%')

%Normalized plots

%Speed
fig.xaclxq_tot_alpha_x=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Arrival Rate')


%Quality
fig.dqaclxq_tot_alpha_x=figure; hold all
%sb=scatter(eq.lx/par.L_I,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Quality')


%Growth
fig.gaclxq_tot_alpha_x=figure; hold all
%sb=scatter(eq.lx/par.L_I,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
xlabel('Proportion To Speed')
ylabel('Growth \%')

% Cell with colors
pcolorc={'k',par.opt.light_blue,par.opt.maroon,par.opt.green};
stylec={':','-','--','-.'};

lgds=[]; lgdq=[]; lgdg=[];
lgdsa=[]; lgdqa=[]; lgdga=[];
lgdstr=cell(1, npar);

%Calibration value
xcal=eq.x;
dqcal=eq.dq;
gcal=eq.g;


for ipar=1:npar

    %Update parameter
    par.alpha_x=alpha_xv(ipar);

    %Changing the speed and quality of labor (Relative to total labor)

    %Arrival rate
    xg=par.chi*lxg.^par.alpha_x;
    x0g=par0.chi*lx0g.^par0.alpha_x;

    %Quality
    dqg=par.lambda*lqg.^par.alpha_q;
    dq0g=par0.lambda*lq0g.^par0.alpha_q;

    %Entry
    xe=par.chi_e*eq.le.^par.alpha_e;
    dqe=par.lambda_e;

    xe0=par0.chi_e*eq0.le.^par0.alpha_e;
    dqe0=par0.lambda_e;

    %Growth
    gg=xe.*dqe+xg.*dqg;
    g0g=xe0.*dqe0+x0g.*dq0g;


    %Proportion relative to total

    %-----------------------------
    %Reallocation of Labor
    %-----------------------------

    %Speed
    figure(fig.xclxq_tot_alpha_x);
    ps(ipar)=plot(lxg/par.L_I,xg,stylec{ipar}); ps(ipar).Color=pcolorc{ipar};

    %Quality
    figure(fig.dqclxq_tot_alpha_x)
    pq(ipar)=plot(lxg/par.L_I,dqg,stylec{ipar}); pq(ipar).Color=pcolorc{ipar};


    %Growth
    figure(fig.gclxq_tot_alpha_x)
    pg(ipar)=plot(lxg/par.L_I,100*gg,stylec{ipar}); pg(ipar).Color=pcolorc{ipar};

    %Legends
    lgdstr{ipar}=['$\alpha_x$=' num2str(par.alpha_x,2)  ];
    lgds=[lgds,ps(ipar)]; lgdq=[lgdq,pq(ipar)];  lgdg=[lgdg,pg(ipar)];

    %-----------------------------
    %Normalizing to calibration
    %-----------------------------

    %Additive normalization
    xag=xg+(xcal-xg(ind));
    dqag=dqg+(dqcal-dqg(ind));
    gag=gg+(gcal-gg(ind));

     %Find max
    [~,indmax]=max(gag);
    plot_max=1;


    %Speed
    figure(fig.xaclxq_tot_alpha_x);
    psa(ipar)=plot(lxg/par.L_I,xag,stylec{ipar}); psa(ipar).Color=pcolorc{ipar};

    %Quality
    figure(fig.dqaclxq_tot_alpha_x)
    pqa(ipar)=plot(lxg/par.L_I,dqag,stylec{ipar}); pqa(ipar).Color=pcolorc{ipar};


    %Growth
    figure(fig.gaclxq_tot_alpha_x)
    pga(ipar)=plot(lxg/par.L_I,100*gag,stylec{ipar}); pga(ipar).Color=pcolorc{ipar};

    if plot_max==1
        plxg=lxg/par.L_I;
        smax=scatter(plxg(indmax),100*gag(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=pcolorc{ipar}; smax.MarkerEdgeColor=pcolorc{ipar};
    end

    %Legends
    if ipar==1
        lgdastr{ipar}=['Calibration $\alpha_x$=' num2str(par.alpha_x,2)  ];
        lgdsa=[lgdsa,psa(ipar)]; lgdqa=[lgdqa,pqa(ipar)];  lgdga=[lgdga,pga(ipar)];
    else
        lgdastr{ipar}=['$\alpha_x$=' num2str(par.alpha_x,2)  ];
        lgdsa=[lgdsa,psa(ipar)]; lgdqa=[lgdqa,pqa(ipar)];  lgdga=[lgdga,pga(ipar)];
    end
end

%Calibration lines and legend

figure(fig.xclxq_tot_alpha_x);
yLimits = get(gca,'YLim');
cx=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgds,lgdstr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="northwest";

figure(fig.dqclxq_tot_alpha_x)
yLimits = get(gca,'YLim');
cq=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdq,lgdstr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

figure(fig.gclxq_tot_alpha_x)
yLimits = get(gca,'YLim');
cg=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdg,lgdstr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

%Normalized
figure(fig.xaclxq_tot_alpha_x);
yLimits = get(gca,'YLim');
cx=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdsa,lgdastr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

figure(fig.dqaclxq_tot_alpha_x)
yLimits = get(gca,'YLim');
cq=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdqa,lgdastr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="best";

figure(fig.gaclxq_tot_alpha_x)
yLimits = get(gca,'YLim');
cg=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,':k');
lgd=legend(lgdga,lgdastr); legend('boxoff'); lgd.EdgeColor='none'; lgd.Location="south";
ylim([0.5,2])



%% Save figures
%par.opt.save_fig=1;
if ~isfield(par,'save_str')
    par.save_str='';
end

if par.opt.save_fig==1
    %alpha_q
    exportgraphics(fig.xclxq_tot_alpha_q,[par.opt.dir_fig 'xclxq_tot_alpha_q' par.save_str '.pdf'])
    exportgraphics(fig.dqclxq_tot_alpha_q,[par.opt.dir_fig 'dqclxq_tot_alpha_q' par.save_str '.pdf'])
    exportgraphics(fig.gclxq_tot_alpha_q,[par.opt.dir_fig 'gclxq_tot_alpha_q' par.save_str '.pdf'])
    exportgraphics(fig.xaclxq_tot_alpha_q,[par.opt.dir_fig 'xaclxq_tot_alpha_q' par.save_str '.pdf'])
    exportgraphics(fig.dqaclxq_tot_alpha_q,[par.opt.dir_fig 'dqaclxq_tot_alpha_q' par.save_str '.pdf'])
    exportgraphics(fig.gaclxq_tot_alpha_q,[par.opt.dir_fig 'gaclxq_tot_alpha_q' par.save_str '.pdf'])

    %alpha_x
    exportgraphics(fig.xclxq_tot_alpha_x,[par.opt.dir_fig 'xclxq_tot_alpha_x' par.save_str '.pdf'])
    exportgraphics(fig.dqclxq_tot_alpha_x,[par.opt.dir_fig 'dqclxq_tot_alpha_x' par.save_str '.pdf'])
    exportgraphics(fig.gclxq_tot_alpha_x,[par.opt.dir_fig 'gclxq_tot_alpha_x' par.save_str '.pdf'])
    exportgraphics(fig.xaclxq_tot_alpha_x,[par.opt.dir_fig 'xaclxq_tot_alpha_x' par.save_str '.pdf'])
    exportgraphics(fig.dqaclxq_tot_alpha_x,[par.opt.dir_fig 'dqaclxq_tot_alpha_x' par.save_str '.pdf'])
    exportgraphics(fig.gaclxq_tot_alpha_x,[par.opt.dir_fig 'gaclxq_tot_alpha_x' par.save_str '.pdf'])

end



end