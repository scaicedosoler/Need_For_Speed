%Untargeted moments test
function untarget_mom_intq_B(par0,eq0,par1,eq1)

par=par0;

%Compute log differences

%Arrival rate
dx0=log(eq0.x_top)-log(eq0.x_bot);
dx1=log(eq1.x_top)-log(eq1.x_bot);
Dx1_0=dx1-dx0;

%Speed: arrival per inventor
dx_l0=log(eq0.x_top/eq0.l_top)-log(eq0.x_bot/eq0.l_bot);
dx_l1=log(eq1.x_top/eq1.l_top)-log(eq1.x_bot/eq1.l_bot);
Dxl1_0=dx_l1-dx_l0;

%Quality
dQ0=log(eq0.Q_top)-log(eq0.Q_bot);
dQ1=log(eq1.Q_top)-log(eq1.Q_bot);
DQ1_0=dQ1-dQ0;

%Value
dV0=log(eq0.V_top)-log(eq0.V_bot);
dV1=log(eq1.V_top)-log(eq1.V_bot);
DV1_0=dV1-dV0;


%% Plots

%Number of periods
tg=[-3,-2,-1,0,1,2,3];
ng=length(tg);


%Arrival Rate
fig_dx=figure; hold all;
p0=plot(tg,dx0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
p1=plot(tg,dx1*(tg>=0),'--D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
plot(tg,zeros(1,ng),'k:')
lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Periods'); ylabel('Log Arrival Rate');

%Speed: arrival per inventor
fig_dx_l=figure; hold all;
p0=plot(tg,dx_l0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
p1=plot(tg,dx_l1*(tg>=0),'--D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
plot(tg,zeros(1,ng),'k:')
lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Periods'); ylabel('Log Speed: Arrival Rate per Inventor');
legend('off')

%Quality
fig_dQ=figure; hold all;
p0=plot(tg,dQ0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
p1=plot(tg,dQ1*(tg>=0),'--D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
plot(tg,zeros(1,ng),'k:')
lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Periods'); ylabel('Log Quality');
legend('off')

%Value
fig_dV=figure; hold all;
p0=plot(tg,dV0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
p1=plot(tg,dV1*(tg>=0),'--D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
plot(tg,zeros(1,ng),'k:')
lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
xlabel('Periods'); ylabel('Log Value');
legend('off')


%% Save figures
%par.opt.save_fig=1;
if ~isfield(par,'save_str')
    par.save_str='';
end

if par.opt.save_fig==1
    exportgraphics(fig_dx,[par.opt.dir_fig 'dx_untarget_mom' par.save_str '.pdf'])
    exportgraphics(fig_dx_l,[par.opt.dir_fig 'dx_l_untarget_mom' par.save_str '.pdf'])
    exportgraphics(fig_dQ,[par.opt.dir_fig 'dQ_untarget_mom' par.save_str '.pdf'])
    exportgraphics(fig_dV,[par.opt.dir_fig 'dV_untarget_mom' par.save_str '.pdf'])

end

%% Data

par.use_data=1;

if ~isfield(par,'use_data')
    par.use_data=0;
end

if par.use_data==1

    %Choose variables
    dx_str='lnnpat_f';
    dx_l_str='lnspeed';
    dQ_str='lf_cit3';
    dV_str='logxi';

    %Choose estimation
    est_str='base';

    %Load data
    load('data_un_mom_c1_none_balanced_gr1090_f.mat')
    %load('data_un_mom_move_no_kmatch_1980_bottom_top_vs_bottom__c1_none_balanced_gr1090_f.mat')
    %par.save_str='_move_no_kmatch';

    ind_t0g=find(data.(dx_str).base.event_year>=-3&data.(dx_str).base.event_year<=3&data.(dx_str).base.pre1985==1);
    ind_t1g=find(data.(dx_str).base.event_year>=-3&data.(dx_str).base.event_year<=3&data.(dx_str).base.pre1985==0);

    %Arrival Rate
    fig_dx_data=figure; hold all;
    fill([tg, fliplr(tg)], [data.(dx_str).base.max95(ind_t0g)', fliplr(data.(dx_str).base.min95(ind_t0g)')], par0.opt.light_blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([tg, fliplr(tg)], [data.(dx_str).base.max95(ind_t1g)', fliplr(data.(dx_str).base.min95(ind_t1g)')], par0.opt.maroon, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    pd0=plot(tg,data.(dx_str).base.estimate(ind_t0g),'--'); pd0.Color=par0.opt.light_blue; pd0.MarkerFaceColor='w'; pd0.MarkerEdgeColor=par0.opt.light_blue;
    pd1=plot(tg,data.(dx_str).base.estimate(ind_t1g),'--'); pd1.Color=par0.opt.maroon; pd1.MarkerFaceColor='w'; pd1.MarkerEdgeColor=par0.opt.maroon;
    p0=plot(tg,dx0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
    p1=plot(tg,dx1*(tg>=0),'-D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
    plot(tg,zeros(1,ng),'k:')
    lgd=legend([pd0,pd1,p0,p1],{['Data ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Data ' num2str(par1.iyear) '-' num2str(par1.iyear+5)],['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Periods'); ylabel('Log Arrival Rate');

    %Speed: arrival per inventor
    fig_dx_l_data=figure; hold all;
    fill([tg, fliplr(tg)], [data.(dx_l_str).base.max95(ind_t0g)', fliplr(data.(dx_l_str).base.min95(ind_t0g)')], par0.opt.light_blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([tg, fliplr(tg)], [data.(dx_l_str).base.max95(ind_t1g)', fliplr(data.(dx_l_str).base.min95(ind_t1g)')], par0.opt.maroon, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    pd0=plot(tg,data.(dx_l_str).base.estimate(ind_t0g),'--'); pd0.Color=par0.opt.light_blue; pd0.MarkerFaceColor='w'; pd0.MarkerEdgeColor=par0.opt.light_blue;
    pd1=plot(tg,data.(dx_l_str).base.estimate(ind_t1g),'--'); pd1.Color=par0.opt.maroon; pd1.MarkerFaceColor='w'; pd1.MarkerEdgeColor=par0.opt.maroon;
    p0=plot(tg,dx_l0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
    p1=plot(tg,dx_l1*(tg>=0),'-D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
    plot(tg,zeros(1,ng),'k:')
    lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Periods'); ylabel('Log Speed: Arrival Rate per Inventor');
    legend('off')

    %Quality
    fig_dQ_data=figure; hold all;
    fill([tg, fliplr(tg)], [data.(dQ_str).base.max95(ind_t0g)', fliplr(data.(dQ_str).base.min95(ind_t0g)')], par0.opt.light_blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([tg, fliplr(tg)], [data.(dQ_str).base.max95(ind_t1g)', fliplr(data.(dQ_str).base.min95(ind_t1g)')], par0.opt.maroon, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    pd0=plot(tg,data.(dQ_str).base.estimate(ind_t0g),'--'); pd0.Color=par0.opt.light_blue; pd0.MarkerFaceColor='w'; pd0.MarkerEdgeColor=par0.opt.light_blue;
    pd1=plot(tg,data.(dQ_str).base.estimate(ind_t1g),'--'); pd1.Color=par0.opt.maroon; pd1.MarkerFaceColor='w'; pd1.MarkerEdgeColor=par0.opt.maroon;
    p0=plot(tg,dQ0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
    p1=plot(tg,dQ1*(tg>=0),'-D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
    plot(tg,zeros(1,ng),'k:')
    lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Periods'); ylabel('Log Quality');
    legend('off')

    %Value
    fig_dV_data=figure; hold all;
    fill([tg, fliplr(tg)], [data.(dV_str).base.max95(ind_t0g)', fliplr(data.(dV_str).base.min95(ind_t0g)')], par0.opt.light_blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([tg, fliplr(tg)], [data.(dV_str).base.max95(ind_t1g)', fliplr(data.(dV_str).base.min95(ind_t1g)')], par0.opt.maroon, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    pd0=plot(tg,data.(dV_str).base.estimate(ind_t0g),'--'); pd0.Color=par0.opt.light_blue; pd0.MarkerFaceColor='w'; pd0.MarkerEdgeColor=par0.opt.light_blue;
    pd1=plot(tg,data.(dV_str).base.estimate(ind_t1g),'--'); pd1.Color=par0.opt.maroon; pd1.MarkerFaceColor='w'; pd1.MarkerEdgeColor=par0.opt.maroon;
    p0=plot(tg,dV0*(tg>=0),'-O'); p0.Color=par0.opt.light_blue; p0.MarkerFaceColor='w'; p0.MarkerEdgeColor=par0.opt.light_blue;
    p1=plot(tg,dV1*(tg>=0),'-D'); p1.Color=par0.opt.maroon; p1.MarkerFaceColor='w'; p1.MarkerEdgeColor=par0.opt.maroon;
    plot(tg,zeros(1,ng),'k:')
    lgd=legend([p0,p1],{['Calibration ' num2str(par0.iyear) '-' num2str(par0.iyear+5)],['Calibration ' num2str(par1.iyear) '-' num2str(par1.iyear+5)]}); lgd.Location='northwest';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Periods'); ylabel('Log Value');
    legend('off')


    %% Triple diff

    % From regressions in movers_out_of_sample_draft8.do
    dx_d01=.485327;
    dx_min95_01=.3385149;
    dx_max95_01=.6321391;

    dxl_d01=-.1908427;
    dxl_min95_01=-.3550386;
    dxl_max95_01=-.0266468;

    dQ_d01=-.205482;
    dQ_min95_01=-.2805382;
    dQ_max95_01=-.1304258;

    dV_d01=-.0135524;
    dV_min95_01=-.3219797;
    dV_max95_01=.2948749;


    %Triple Diff
    %All
    fig_TD_all=figure; hold all;
    px=plot([1,1], [dx_min95_01, dx_max95_01]); px.Color=par.opt.light_blue;
    pxl=plot([2,2], [dxl_min95_01, dxl_max95_01]); pxl.Color=par.opt.light_blue;
    pQ=plot([3,3], [dQ_min95_01, dQ_max95_01]); pQ.Color=par.opt.light_blue;
    pV=plot([4,4], [dV_min95_01, dV_max95_01]); pV.Color=par.opt.light_blue;
    sd=scatter(1:4, [dx_d01, dxl_d01, dQ_d01,dV_d01]); sd.Marker="o"; sd.MarkerFaceColor=par.opt.light_blue; sd.MarkerEdgeColor=par.opt.light_blue;
    s=scatter(1:4, [Dx1_0, Dxl1_0, DQ1_0,DV1_0]); s.Marker="diamond"; s.MarkerFaceColor="w"; s.MarkerEdgeColor=par.opt.maroon;
    plot([0,5], [0,0],":k")
    xticks(1:4); xticklabels({'Arrival Rate','Speed', 'Quality','Value'})
    ylabel('Move $\times$  (Small $\to$ Large) $\times ~I_{2010-2015}$')
    lgd=legend([s,sd], {'Model','Data'}); lgd.Box="off";

    %Arrival rate, quality and value
    fig_TD_x_Q_V=figure; hold all;
    px=plot([1,1], [dx_min95_01, dx_max95_01]); px.Color=par.opt.light_blue;
    pQ=plot([2,2], [dQ_min95_01, dQ_max95_01]); pQ.Color=par.opt.light_blue;
    pV=plot([3,3], [dV_min95_01, dV_max95_01]); pV.Color=par.opt.light_blue;
    sd=scatter(1:3, [dx_d01, dQ_d01,dV_d01]); sd.Marker="o"; sd.MarkerFaceColor=par.opt.light_blue; sd.MarkerEdgeColor=par.opt.light_blue;
    s=scatter(1:3, [Dx1_0, DQ1_0,DV1_0]); s.Marker="diamond"; s.MarkerFaceColor="w"; s.MarkerEdgeColor=par.opt.maroon;
    plot([0,4], [0,0],":k")
    xticks(1:4); xticklabels({'Arrival Rate','Quality','Value'})
    ylabel('Move $\times$ (Small $\to$ Large) $\times ~I_{2010-2015}$')
    lgd=legend([s,sd], {'Model','Data'}); lgd.Box="off";

    %Arrival rate and quality
    fig_TD_x_Q=figure; hold all;
    px=plot([1,1], [dx_min95_01, dx_max95_01]); px.Color=par.opt.light_blue;
    pQ=plot([2,2], [dQ_min95_01, dQ_max95_01]); pQ.Color=par.opt.light_blue;
    sd=scatter(1:2, [dx_d01, dQ_d01]); sd.Marker="o"; sd.MarkerFaceColor=par.opt.light_blue; sd.MarkerEdgeColor=par.opt.light_blue;
    s=scatter(1:2, [Dx1_0, DQ1_0]); s.Marker="diamond"; s.MarkerFaceColor="w"; s.MarkerEdgeColor=par.opt.maroon;
    plot([0,3], [0,0],":k")
    xticks(1:4); xticklabels({'Arrival Rate','Quality'})
    ylabel('Move $\times$ (Small $\to$ Large)  $\times ~I_{2010-2015}$')
    lgd=legend([s,sd], {'Model','Data'}); lgd.Box="off";

    if par.opt.save_fig==1
        exportgraphics(fig_TD_all,[par.opt.dir_fig 'TD_all_untarget_mom' par.save_str '.pdf'])
        exportgraphics(fig_TD_x_Q_V,[par.opt.dir_fig 'TD_x_Q_V_untarget_mom' par.save_str '.pdf'])
        exportgraphics(fig_TD_x_Q,[par.opt.dir_fig 'TD_x_Q_untarget_mom' par.save_str '.pdf'])
    end

end



end





