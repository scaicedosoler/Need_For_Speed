%Untargeted moments test
function untarget_mom_intq_B(par0,eq0,par1,eq1)

par=par0;

%Compute log differences

%Arrival rate
dx0=log(eq0.x_top)-log(eq0.x_bot);
dx1=log(eq1.x_top)-log(eq1.x_bot);

%Speed: arrival per inventor
dx_l0=log(eq0.x_top/eq0.l_top)-log(eq0.x_bot/eq0.l_bot);
dx_l1=log(eq1.x_top/eq1.l_top)-log(eq1.x_bot/eq1.l_bot);

%Quality
dQ0=log(eq0.Q_top)-log(eq0.Q_bot);
dQ1=log(eq1.Q_top)-log(eq1.Q_bot);

%Value
dV0=log(eq0.V_top)-log(eq0.V_bot);
dV1=log(eq1.V_top)-log(eq1.V_bot);

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

    ind_t0g=find(data.(dx_str).base.event_year>=-3&data.(dx_str).base.event_year<=3&data.(dx_str).base.pre1985==1);
    ind_t1g=find(data.(dx_str).base.event_year>=-3&data.(dx_str).base.event_year<=3&data.(dx_str).base.pre1985==0);

    %Choose version of estimation

    


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

    if par.opt.save_fig==1
    exportgraphics(fig_dx_data,[par.opt.dir_fig 'dx_data_untarget_mom' par.save_str '.pdf'])
    exportgraphics(fig_dx_l_data,[par.opt.dir_fig 'dx_l_data_untarget_mom' par.save_str '.pdf'])
    exportgraphics(fig_dQ_data,[par.opt.dir_fig 'dQ_data_untarget_mom' par.save_str '.pdf'])
    exportgraphics(fig_dV_data,[par.opt.dir_fig 'dV_data_untarget_mom' par.save_str '.pdf'])
    end

end





end





