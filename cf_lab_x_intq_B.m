%Counterfactual: Reallocating labor with external effects
function cf_lab_x_intq_B(par,eq,par0,eq0)



%-------------------------------------------
% 1. Changing the speed and quality of labor--With external effects
%-------------------------------------------

%External effect elasticity
zetav=[0, 0.12,0.25,0.5]; 
nv=length(zetav);

ng=1000;

%Incumbent labor
linc=eq.lx+eq.lq;
linc0=eq0.lx+eq0.lq;

lxg=linspace(0.01*linc,0.99*linc,ng);
lx0g=linspace(0.01*linc0,0.99*linc0,ng);

plxg=lxg./linc;
plx0g=lx0g./linc0;

lqg=(1-plxg)*linc;
lq0g=(1-plx0g)*linc0;

%Arrival rate
xg=par.chi*lxg.^par.alpha_x;
x0g=par0.chi*lx0g.^par0.alpha_x;

%Quality
dqg=par.lambda.*lqg.^par.alpha_q;
dq0g=par0.lambda*lq0g.^par0.alpha_q;

%Entry
xe=par.chi_e*eq.le.^par.alpha_e;
dqe=par.lambda_e;

xe0=par0.chi_e*eq0.le.^par0.alpha_e;
dqe0=par0.lambda_e;

%Growth
gg=xe.*dqe+xg.*dqg;


% Cell with colors
pcolorc={'k',par.opt.light_blue,par.opt.maroon,par.opt.green};
stylec={':','-','--','-.'};


% Initialize plots

%Growth (initial period)
fig.g_x0clxq=figure; hold all ;
xlabel('Proportion of Labor To Speed')
ylabel('Growth \%')

%Growth (final period)
fig.g_xcmaxlxq=figure; hold all ;
xlabel('Proportion of Labor To Speed')
ylabel('Growth \%')

%Growth (initial period using final speed)
fig.g_x0xclxq=figure; hold all ;
xlabel('Proportion of Labor To Speed')
ylabel('Growth \%')


lgd_str=cell(1, nv);
lgd_g=[]; lgd_g0=[]; lgd_g0x=[];


%% With external effect

for i=1:nv

    %Update zeta
    zeta=zetav(i);

    %Quality
    dq_x0g=par0.lambda*(x0g./eq0.x).^(-zeta).*lq0g.^par0.alpha_q;
    dq_xg=par.lambda*(xg/eq.x).^(-zeta).*lqg.^par.alpha_q;
    dq_x0xg=par0.lambda*(xg./eq0.x).^(-zeta).*lq0g.^par0.alpha_q;
    
    %Growth
    g_xg=xe.*dqe+xg.*dq_xg;
    g_x0g=xe0.*dqe0+x0g.*dq_x0g;
    g_x0xg=xe0.*dqe0+x0g.*dq_x0xg;

    %Compute the maximum
    [~,indmax_x]=max(g_xg);

    %Plot
    %Growth (initial period)
    figure(fig.g_x0clxq)
    %p0(i)=plot(plx0g,100*g_x0g,stylec{i}); p0(i).Color=pcolorc{i};
    p0(i)=plot(lx0g/par0.L_I,100*g_x0g,stylec{i}); p0(i).Color=pcolorc{i};
    %smax0=scatter(plx0g(indmax_x0),100*g_x0g(indmax_x0),80);  smax0.Marker='pentagram'; smax0.MarkerFaceColor=pcolorc{i}; smax0.MarkerEdgeColor=pcolorc{i};

    %Growth (final period)
    figure(fig.g_xcmaxlxq)
    %p1(i)=plot(plxg,100*g_xg,stylec{i}); p1(i).Color=pcolorc{i};
    p1(i)=plot(lxg/par.L_I,100*g_xg,stylec{i}); p1(i).Color=pcolorc{i};    
    smax=scatter(lxg(indmax_x)/par.L_I,100*g_xg(indmax_x),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=pcolorc{i}; smax.MarkerEdgeColor=pcolorc{i};

    %Growth (initial period using final speed)
    figure(fig.g_x0xclxq)
    %p0(i)=plot(plx0g,100*g_x0g,stylec{i}); p0(i).Color=pcolorc{i};
    p0x(i)=plot(lx0g/par0.L_I,100*g_x0xg,stylec{i}); p0x(i).Color=pcolorc{i};

    %Legends
    if i==1
        lgd_str{i}=['Baseline $\zeta$=' num2str(zeta,2)  ];
    end
   
    if i==2
        lgd_str{i}=['Calibration $\zeta$=' num2str(zeta,2)  ];
    end

    if i>2
        lgd_str{i}=['$\zeta$=' num2str(zeta,2)  ];
    end

    lgd_g0=[lgd_g0,p0(i)];
    lgd_g=[lgd_g,p1(i)]; 
    lgd_g0x=[lgd_g0x,p0x(i)];

end


%Calibration lines and legend

%Growth (initial period)
figure(fig.g_x0clxq);
%s0=scatter(eq0.lx/linc0,100*eq0.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor='k';
s0=scatter(eq0.lx/par0.L_I,100*eq0.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor='k';
ylim([1.5,1.68]); xlim([0.5,1])
yLimits = get(gca,'YLim');
cx=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,'k');
lgd=legend(lgd_g0,lgd_str); legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; lgd.NumColumns=1; 

%Growth (final period)
figure(fig.g_xcmaxlxq);
%s0=scatter(eq.lx/linc,100*eq.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor='k';
s0=scatter(eq.lx/par.L_I,100*eq.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor='k';
ylim([1.2,1.4]); xlim([0.1,1])
lgd=legend(lgd_g,lgd_str); legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="northwest"; lgd.NumColumns=1; 

%Growth (initial period using final speed)
figure(fig.g_x0xclxq);
%s0=scatter(eq.lx/linc,100*eq.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor='k';
s0=scatter(eq0.lx/par0.L_I,100*eq0.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor='k';
ylim([1.2,1.68]); xlim([0.5,1])
yLimits = get(gca,'YLim');
cx=plot([eq.lx/par.L_I,eq.lx/par.L_I],yLimits,'k');
lgd=legend(lgd_g0x,lgd_str); legend('boxoff'); lgd.EdgeColor='none';  lgd.Location="best"; lgd.NumColumns=1; 


%Save figures
%par.opt.save_fig=1;
if ~isfield(par,'save_str')
    par.save_str='';
end

if par.opt.save_fig==1

    exportgraphics(fig.g_x0clxq,[par.opt.dir_fig 'g_x0clxq_tot_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.g_xcmaxlxq,[par.opt.dir_fig 'g_xcmaxlxq_tot_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.g_x0xclxq,[par.opt.dir_fig 'g_x0xclxq_tot_intq_B' par.save_str '.pdf'])
    
end

 

end