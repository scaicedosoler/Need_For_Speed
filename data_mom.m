% Need for Speed
% Caicedo-Pearce
% Moments from data
% Spring 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all;
% clc, clear all;
%close all;

%Figure format
par=fig_format();

%Path
addpath('./mat/')
par.opt.dir_mat='./mat/';
par.opt.dir_fig='./figures/';
par.opt.dir_tab='./tables/';

%Saving options
par.opt.save_fig=0;
par.opt.save_tab=0;

%Running options
plot_data=0;
save_data=1;
save_fig=0;
save_str='_ipc3'; %'_ipc3_ns'; % '_ns'; %  '_im'; %'_ns_im'; %

warning('off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upload data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir= '/Users/nicolasfajardo/Dropbox/Caicedo_Pearce/RA - Nicolas/tasks/Task 28 - Model and Matlab Code/Quality_Quantity/data/'; %'..\data\'; %'D:\Dropbox\santiago\Research\Caicedo_Pearce\RA - Nicolas\tasks\Task 15 - Moments for Exogenous Types with Multiple Skills\'; %'..\data\'; 

filename =[dir 'mom_data' save_str '.xlsx']; % [dir 'mom_data_SC.xlsx']; %[dir 'data_mom_example.xlsx'];

%Create data structure
data = readtable(filename);

%Additional variable definitions
data.pat_90=data.pat_tot-data.pat_10;
data.pat_int_90=data.pat_int-data.pat_int_10;
data.pat_ext_90=data.pat_ext-data.pat_ext_10;
data.pat_inc_90=data.pat_inc-data.pat_inc_10;

data.inv_90=data.inv_tot-data.inv_10;
data.inv_inc_90=data.inv_inc-data.inv_inc_10;
data.inv_int_90=data.inv_int-data.inv_int_10;
data.inv_ext_90=data.inv_ext-data.inv_ext_10;

data.inv_int_hs_90=data.inv_int_hs-data.inv_int_hs_10;
data.inv_ext_hs_90=data.inv_ext_hs-data.inv_ext_hs_10;

%----------------------------------------
%Create structure with moments to match
%----------------------------------------
dy=4;    %Number of years for average

%Base year for comparison moments
indb=find(data.year>=1980 & data.year<=1980+dy);

for iy=1978:dy:(2014-dy)

    ind=find(data.year>=iy & data.year<=iy+dy);

    dmom.iyear=iy;
    dmom.fyear=iy+dy;

    %-----------
    % Patents
    %-----------

    %Concentration (patents)
    dmom.pat_10=mean(data.pat_10(ind)./data.pat_tot(ind));

    %Share of incumbent innovations
    dmom.pat_inc=mean(data.pat_inc(ind)./data.pat_tot(ind));

    %Share of entrants innovations
    dmom.pat_ent=mean(data.pat_ent(ind)./data.pat_tot(ind));

    %Share of internal innovations
    dmom.pat_int=mean(data.pat_int(ind)./data.pat_tot(ind));

    %Share of external innovations
    dmom.pat_ext=mean(data.pat_ext(ind)./data.pat_tot(ind));

    %Share of internal innovations in top 10 patents
    dmom.pat_int_10=mean(data.pat_int_10(ind)./data.pat_10(ind));

    %Share of internal innovations in bottom 90 patents
    dmom.pat_int_90=mean(data.pat_int_90(ind)./data.pat_90(ind));

    %Share of top 10 patents in total internal innovations
    dmom.pat_10_int=mean(data.pat_int_10(ind)./data.pat_int(ind));

    %Share of top 10 patents in total external innovations
    dmom.pat_10_ext=mean(data.pat_ext_10(ind)./data.pat_ext(ind));

    %Patents per firm
    dmom.pat_f=mean(data.pat_tot./data.f_tot);
    dmom.pat_f_10=mean(data.pat_10./data.f_10);
    dmom.pat_f_90=mean(data.pat_90./data.f_90);

    %Patents relative to base year
    dmom.pat_rb=sum(data.pat_tot(ind))/sum(data.pat_tot(indb));

    %Patents per firm relative to base year
    dmom.pat_f_rb=sum(data.pat_tot(ind)./data.f_tot(ind))./sum(data.pat_tot(indb)./data.f_tot(indb));

    dmom.pat_f_10_rb=sum(data.pat_10(ind)./data.f_10(ind))./sum(data.pat_10(indb)./data.f_10(indb));
    dmom.pat_f_90_rb=sum(data.pat_90(ind)./data.f_90(ind))./sum(data.pat_90(indb)./data.f_90(indb));

    %-----------
    % Inventors
    %-----------

    %Concentration (inventors)
    dmom.inv_10=mean(data.inv_10(ind)./data.inv_tot(ind));

    %Share of inventors in incumbent innovations
    dmom.inv_inc=mean(data.inv_inc(ind)./data.inv_tot(ind));

    %Share of inventors in entrant innovations
    dmom.inv_ent=mean(data.inv_ent(ind)./data.inv_tot(ind));

    %Share of inventors in internal innovations
    dmom.inv_int=mean(data.inv_int(ind)./data.inv_tot(ind));

    %Share of external innovations
    dmom.inv_ext=mean(data.inv_ext(ind)./data.inv_tot(ind));

    %Share of inventors in internal innovations in top 10
    dmom.inv_int_10=mean(data.inv_int_10(ind)./data.inv_10(ind));

    %Share of inventors in internal innovations in bottom 90
    dmom.inv_int_90=mean(data.inv_int_90(ind)./data.inv_90(ind));

    %Share of top 10 inventors in total internal inventors
    dmom.inv_10_int=mean(data.inv_int_10(ind)./data.inv_int(ind));

    %Share of top 10 inventors in total external innovations
    dmom.inv_10_ext=mean(data.inv_ext_10(ind)./data.inv_ext(ind));

    %Share of high skill inventors in top 10
    dmom.inv_hs_10=mean(data.inv_hs_10(ind)./data.inv_tot_hs(ind));

    %Share of incumbents inventors that are high-skill
    dmom.inv_inc_hs=mean(data.inv_inc_hs(ind)./data.inv_inc(ind));

    %Share of entrant inventors that are high-skill
    dmom.inv_ent_hs=mean(data.inv_ent_hs(ind)./data.inv_ent(ind));

    %Share of internal inventors that are high skill
    dmom.inv_int_hs=mean(data.inv_int_hs(ind)./data.inv_int(ind));

    %Share of external inventors that are high skill
    dmom.inv_ext_hs=mean(data.inv_ext_hs(ind)./data.inv_ext(ind));

    %Share of high-skill that are incumbent inventors
    dmom.inv_hs_inc=mean(data.inv_inc_hs(ind)./data.inv_tot_hs(ind));

    %Share of high-skill that are entrant inventors
    dmom.inv_hs_ent=mean(data.inv_ent_hs(ind)./data.inv_tot_hs(ind));

    %Share of high-skill that are internal inventors
    dmom.inv_hs_int=mean(data.inv_int_hs(ind)./data.inv_tot_hs(ind));

    %Share of high-skill that are external inventors
    dmom.inv_hs_ext=mean(data.inv_ext_hs(ind)./data.inv_tot_hs(ind));

    %Share of high-skill inventors in internal top 10 patents
    dmom.inv_int_hs_10=mean(data.inv_int_hs_10(ind))/mean(data.inv_int_10(ind));

    %Share of high-skill inventors in external top 10 patents
    dmom.inv_ext_hs_10=mean(data.inv_ext_hs_10(ind))/mean(data.inv_ext_10(ind));

    %Share of high-skill inventors in internal bottom 90 patents
    dmom.inv_int_hs_90=mean(data.inv_int_hs_90(ind))/mean(data.inv_int_90(ind));

    %Share of high-skill inventors in external bottom 90 patents
    dmom.inv_ext_hs_90=mean(data.inv_ext_hs_90(ind))/mean(data.inv_ext_90(ind));

    %Inventors relative to base year
    dmom.inv_rb=sum(data.inv_tot(ind))/sum(data.inv_tot(indb));

    %----------------------
    % Patent quality 
    %----------------------
    %Total
    dmom.q=mean(data.q_tot(ind));

    %Internal
    dmom.q_int=mean(data.q_int(ind));

    %External
    dmom.q_ext=mean(data.q_ext(ind));

    %Incumbents
    dmom.q_inc=mean(data.q_inc(ind));

    %Entrants
    dmom.q_ent=mean(data.q_ent(ind));

    %Internal in top 10
    dmom.q_int_10=mean(data.q_int_10(ind));

    %External in top 10
    dmom.q_ext_10=mean(data.q_ext_10(ind));

    %Internal in bottom 90
    dmom.q_int_90=mean(data.q_int_90(ind));

    %External in bottom 90
    dmom.q_ext_90=mean(data.q_ext_90(ind));

    %Relative patent quality internal/external top 10
    dmom.q_int_ext_10=dmom.q_int_10/dmom.q_ext_10;

    %Relative patent quality internal/external bottom 90
    dmom.q_int_ext_90=dmom.q_int_90/dmom.q_ext_90;

    %Relative patent quality entrant/incumbents
    dmom.q_ent_inc=dmom.q_ent/dmom.q_inc;

    %Without units
    %Internal in top 10 over qbar
    dmom.q_int_10_qbar=mean(data.q_int_10(ind)./data.q_tot(ind));

    %External in top 10
    dmom.q_ext_10_qbar=mean(data.q_ext_10(ind)./data.q_tot(ind));

    %Internal in bottom 90
    dmom.q_int_90_qbar=mean(data.q_int_90(ind)./data.q_tot(ind));

    %External in bottom 90
    dmom.q_ext_90_qbar=mean(data.q_ext_90(ind)./data.q_tot(ind));

    %Relative quality top 10
    dmom.q_10_90=mean(data.q_10(ind)./data.q_90(ind));

    %----------------------------------------
    % Ratios: patents, quality and value per inventor
    %----------------------------------------
    
    %Patents per inventor
    dmom.pat_inv=mean(data.pat_tot(ind)./data.inv_tot(ind));

    %Inventors per patent
    dmom.inv_pat=mean(data.inv_tot(ind)./data.pat_tot(ind));

    %Arrival rate over unique inventors
    dmom.x_l=mean(data.pat_tot(ind)./data.inv_uniq(ind));

    %Quality per inventor
    dmom.lf_cit3_inv=mean(data.logf_cit3(ind)./data.inv_tot(ind));
    dmom.lf_cit3_res_inv=mean(data.logf_cit3_res_ipc_ts(ind)./data.inv_tot(ind));

    %Value per inventor
    dmom.lxi_inv=mean(data.logxi_real(ind)./data.inv_tot(ind)); 
    dmom.lxi_res_inv=mean(data.logxi_real_res_ipc_ts(ind)./data.inv_tot(ind));

    %Value per citation
    dmom.lxi_lf_cit3=mean(data.logxi_real(ind)./data.logf_cit3(ind));
    dmom.lxi_lf_cit3_res=mean(data.logxi_real_res_ipc_ts(ind)./data.logf_cit3_res_ipc_ts(ind));
    dmom.xi_lf_cit3=mean(data.xi_real(ind)./data.logf_cit3(ind));
    dmom.xi_lf_cit3_res=mean(data.xi_real(ind)./data.logf_cit3_res_ipc_ts(ind));

    %----------------------
    % Inventor transitions 
    %----------------------

    %Probability low skill --> high skill
    dmom.pr_ls_hs=mean(data.inv_ls_np_hs(ind)./data.inv_ls_np(ind));

    %Probability high skill --> low skill 
    dmom.pr_hs_ls=mean(data.inv_hs_np_ls(ind)./data.inv_hs_np(ind));

    %Probability low skill --> high skill in top 10
    dmom.pr_ls_hs_10=mean(data.inv_ls_np_hs_10(ind)./data.inv_ls_np_10(ind));

    %Probability high skill --> low skill in top 10
    dmom.pr_hs_ls_10=mean(data.inv_hs_np_ls_10(ind)./data.inv_hs_np_10(ind));

    %Probability low skill --> high skill in bottom 90
    dmom.pr_ls_hs_90=mean(data.inv_ls_np_hs_90(ind)./data.inv_ls_np_90(ind));

    %Probability high skill --> low skill in bottom 90
    dmom.pr_hs_ls_90=mean(data.inv_hs_np_ls_90(ind)./data.inv_hs_np_90(ind));

    %Probability low skill --> high skill in entrants
    dmom.pr_ls_hs_ent=mean(data.inv_ls_np_hs_ent(ind)./data.inv_ls_np_ent(ind));

    %Probability high skill --> low skill in entrants
    dmom.pr_hs_ls_ent=mean(data.inv_hs_np_ls_ent(ind)./data.inv_hs_np_ent(ind));

    %Computed only for final year

    indf=find(data.year>=iy-dy & data.year<=iy); %For death rate take 5 years prior

    %Death rate of high skill
    dmom.pr_hs_out=mean(data.inv_hs_lp(indf)./data.inv_tot_hs(indf));

    %Death rate of low skill
    dmom.pr_ls_out=mean(data.inv_ls_lp(indf)./data.inv_tot_ls(indf));

    %Death rate of high skill in top 10
    dmom.pr_hs_out_10=mean(data.inv_hs_lp_10(indf)./data.inv_hs_10(indf));

    %Death rate of high skill in bottom 90
    dmom.pr_hs_out_90=mean(data.inv_hs_lp_90(indf)./data.inv_hs_90(indf));

    %Death rate of high skill in entrant
    dmom.pr_hs_out_ent=mean(data.inv_hs_lp_ent(indf)./data.inv_ent_hs(indf));

    %Death rate of low skill in top 10
    dmom.pr_ls_out_10=mean(data.inv_ls_lp_10(indf)./data.inv_ls_10(indf));

    %Death rate of low skill in bottom 90
    dmom.pr_ls_out_90=mean(data.inv_ls_lp_90(indf)./data.inv_ls_90(indf));

    %Death rate of low skill in entrant
    dmom.pr_ls_out_ent=mean(data.inv_ls_lp_ent(indf)./data.inv_ent_ls(indf));

    %----------------------
    % Firm transitions 
    %----------------------

    %Probability bottom 90 --> top 10
    dmom.pr_90_10=mean(data.firm_90_np_10(ind)./data.firm_90_np(ind));

    %Probability top 10 --> bottom 90
    dmom.pr_10_90=mean(data.firm_10_np_90(ind)./data.firm_10_np(ind));

    %Probability bottom 90 --> out
    dmom.pr_90_out=mean(data.firm_90_lp(ind)./(data.firm_90_np(ind)+data.firm_90_lp(ind))); %Check this definition- May 2022 -Not at the firm level

    %Probability top 10 --> out
    dmom.pr_10_out=mean(data.firm_10_lp(ind)./(data.firm_10_np(ind)+data.firm_10_lp(ind))); %Check this definition- May 2022 -Not at the firm level

    %----------------------
    % Firm outcomes
    %----------------------
    dmom.pat_val_sales=mean(data.pat_val_sales_nominal(ind));
    dmom.pat_val_sales_r=mean(data.pat_val_sales_real(ind));

    dmom.pat_val_sales_macro=mean(data.pat_val_sales_nominal_macro(ind));
    dmom.pat_val_sales_r_macro=mean(data.pat_val_sales_real_macro(ind));

    dmom.pat_val_sales_win=mean(data.pat_val_sales_nominal_win(ind));
    dmom.pat_val_sales_r_win=mean(data.pat_val_sales_real_win(ind));

    %From literature

    % Profitability from Akcigit and Kerr (2018)
    dmom.prof_AK=0.109;
    % % Average growth from Akcigit and Kerr (2018)
    % dmom.g=0.02;% 0.01; %
    % Entry rate from Akcigit and Kerr (2018)
    dmom.xe_AK=0.058;

    %Klenow et al growth rate
    if iy<1994 %ACA 1995
        dmom.g=0.0166; %For 1983-1993
    end

    if iy>1994 %ACA 1995
        dmom.g=0.0132; %For 2003-2013
    end


    %Additional Moments

    if iy<=1994 %ACA 1995

        %Entrant ratio in top10/bottom90
        dmom.ent_10_90=0.535;
        %Average private value
        dmom.xi_10_90=1.41;
        %Private value/ Patent Quality growth
        dmom.r_lxi_lq=1; %Normalization
        dmom.r_lxi_lq_a=1; %Normalization

        % %Arrival date over unique inventors
        % dmom.x_l=0.77; %Approximation from chart...must compute it correctly : Feb 2024
    end

    if iy>1994 %ACA 1995
        %Entrant ratio in top10/bottom90
        dmom.ent_10_90=0.405;
        %Average private value
        dmom.xi_10_90=1.17;
        %Private value/ Patent Quality growth: see data_mom_xi_f_cit3.xlsx
        dmom.r_lxi_lq=1.70; %Option 1: Ratio (10's/80's) log(xi)/log(fcit3)  residualized gipc3, tsize 
        dmom.r_lxi_lq_a=1.51; %Option 2: Ratio (10's/80's) log(xi)/log(fcit3)  

        % %Arrival rate over unique inventors
        % dmom.x_l=0.67; %Approximation from chart...must compute it correctly : Feb 2024

    end

    % Save data

    if save_data==1
        save([par.opt.dir_mat 'data_mom' save_str '_' num2str(iy) '_' num2str(iy+dy) '.mat'],'data', 'dmom')
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_data==1
    
    %Plot options
    iyear=1980;
    fyear=2014;
    iiy=find(data.year==iyear);
    ify=find(data.year==fyear);

    %----------------------------
    % Patents
    %----------------------------
    
    %Concentration of Innovation
    fig_con=figure; hold all;
    p1=plot(data.year(iiy:ify),data.pat_10(iiy:ify)./data.pat_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of Patents by Top 10 \%');
    xlim([iyear,fyear])
    
    %Share of incumbent patents 
    fig_pat_inc=figure; hold all;
    p1=plot(data.year(iiy:ify),data.pat_inc(iiy:ify)./data.pat_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of Patents by Incumbents');
    xlim([iyear,fyear])

    % Share of internal and external patents
    fig_pat_int_ext=figure; hold all;
    p1=plot(data.year(iiy:ify),data.pat_int(iiy:ify)./data.pat_inc(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.pat_ext(iiy:ify)./data.pat_inc(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Patents');
    xlim([iyear,fyear])

    % Share of internal for top 10 and bottom 90
    data.pat_inc_90=data.pat_inc-data.pat_inc_10;
    data.pat_int_90=data.pat_int-data.pat_int_10;

    fig_pat_int_10_90=figure; hold all;
    p1=plot(data.year(iiy:ify),data.pat_int_10(iiy:ify)./data.pat_inc_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.pat_int_90(iiy:ify)./data.pat_inc_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Internal Patents');
    xlim([iyear,fyear])

    % Share of incumbent for top 10
    fig_pat_inc_10=figure; hold all;
    p1=plot(data.year(iiy:ify),data.pat_inc_10(iiy:ify)./data.pat_inc(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of Incumbent Patents by Top 10');
    xlim([iyear,fyear])
    
    %----------------------------
    % Inventors
    %----------------------------
    
    %Concentration of Inventors
    fig_con_inv=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_10(iiy:ify)./data.inv_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of Inventors in Top 10 \%');
    xlim([iyear,fyear])
    
    %Share of inventors in incumbent firms
    fig_inv_inc=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_inc(iiy:ify)./data.inv_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of Inventors in Incumbents');
    xlim([iyear,fyear])
        
    % Share of inventors in internal and external patents
    fig_inv_int_ext=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int(iiy:ify)./data.inv_inc(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext(iiy:ify)./data.inv_inc(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Inventors');
    xlim([iyear,fyear])

    % Share of high skill in top 10 and bottom 90 
    data.inv_90=data.inv_tot-data.inv_10;

    fig_inv_10_90_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_hs_10(iiy:ify)./data.inv_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_hs_90(iiy:ify)./data.inv_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors');
    xlim([iyear,fyear])

    % Share of inventors in internal from top 10 and bottom 90
    data.inv_int_90=data.inv_int-data.inv_int_10;
    data.inv_ext_90=data.inv_ext-data.inv_ext_10;
    
    fig_inv_int_10_90=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_10(iiy:ify)./data.inv_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_int_90(iiy:ify)./data.inv_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Internal Inventors');
    xlim([iyear,fyear])

    % Share of inventors in external from top 10 and bottom 90
    fig_inv_ext_10_90=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_ext_10(iiy:ify)./data.inv_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_90(iiy:ify)./data.inv_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of External Inventors');
    xlim([iyear,fyear])

    % Share of inventors in internal and external from  top 10
    fig_inv_int_ext_10=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_10(iiy:ify)./data.inv_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_10(iiy:ify)./data.inv_10(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Inventors in Top 10 \%');

    % Share of inventors in internal and external from  bottom 90
    fig_inv_int_ext_90=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_90(iiy:ify)./data.inv_90(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_90(iiy:ify)./data.inv_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Inventors in Bottom 90');
    xlim([iyear,fyear])

    %High skill, low-skill and entrant inventors
    fig_inv_tot_hs_ls_einv=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_tot_hs(iiy:ify)./data.inv_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_tot_ls(iiy:ify)./data.inv_tot(iiy:ify),'--'); p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_tot_einv(iiy:ify)./data.inv_tot(iiy:ify),'--'); p3.Color=par.opt.green;
    lgd=legend([p1,p2,p3],{'High-Skill','Low-Skill', 'New'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of Inventors'); %ylim([0.49,0.51])
    xlim([iyear,fyear])

    %Share of high skill inventors in top 10 
    fig_inv_hs_10=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_hs_10(iiy:ify)./data.inv_tot_hs(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Inventors in Top 10 \%');
    xlim([iyear,fyear])
    
    % Share of high skill in incumbents
    fig_inv_inc_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_inc_hs(iiy:ify)./data.inv_inc(iiy:ify)); p1.Color=par.opt.light_blue;
    ylim([round(0.99*min(data.inv_inc_hs(iiy:ify)./data.inv_inc(iiy:ify)),2),round(1.01*max(data.inv_inc_hs(iiy:ify)./data.inv_inc(iiy:ify)),2)])
    xlabel('Year'); ylabel('Share of High Skill Incument Inventors'); 
    xlim([iyear,fyear])
    
    % Share of high skill in incumbents (out of hs)
    fig_inv_hs_inc=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_inc_hs(iiy:ify)./data.inv_tot_hs(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Inventors in Incumbent');
    
    % Share of high skill in entrants
    fig_inv_ent_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_ent_hs(iiy:ify)./data.inv_ent(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Entrant Inventors');
    xlim([iyear,fyear])
    
    % Share of high skill in entrants (out of hs)
    fig_inv_hs_ent=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_ent_hs(iiy:ify)./data.inv_tot_hs(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Inventors in Entrants');
    xlim([iyear,fyear])
    
    % Share of high skill in internal and external
    fig_inv_int_ext_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs(iiy:ify)./data.inv_int(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs(iiy:ify)./data.inv_ext(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors');
    xlim([iyear,fyear])

    % Share of high skill in internal, external, entrants
    fig_inv_int_ext_ent_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs(iiy:ify)./data.inv_int(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs(iiy:ify)./data.inv_ext(iiy:ify),'--'); p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_ent_hs(iiy:ify)./data.inv_ent(iiy:ify),'--'); p3.Color=par.opt.green;
    lgd=legend([p1,p2,p3],{'Internal','External','Entrant'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors');
    xlim([iyear,fyear])

    % Share of high skill in internal 
    fig_inv_int_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs(iiy:ify)./data.inv_int(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Inventors in Internal');
    xlim([iyear,fyear])
    
    % Share of high skill in external
    fig_inv_ext_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_ext_hs(iiy:ify)./data.inv_ext(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Inventors in External');
    xlim([iyear,fyear])
    
    % Share of high skill in internal and external (out of hs)
    fig_inv_hs_int_ext_ent=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs(iiy:ify)./data.inv_tot_hs(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs(iiy:ify)./data.inv_tot_hs(iiy:ify),'--'); p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_ent_hs(iiy:ify)./data.inv_tot_hs(iiy:ify),':'); p3.Color=par.opt.green;
    lgd=legend([p1,p2,p3],{'Internal','External','Entrants'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors');
    xlim([iyear,fyear])
    
    %Share of hs internal in 10 and 90
    data.inv_int_hs_90=data.inv_int_hs-data.inv_int_hs_10;
    data.inv_int_90=data.inv_int-data.inv_int_10;
    
    fig_inv_int_10_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs_10(iiy:ify)./data.inv_int_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_int_hs_90(iiy:ify)./data.inv_int_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Internal Inventors');
    xlim([iyear,fyear])
    
    %Share of hs external in 10 and 90
    data.inv_ext_hs_90=data.inv_ext_hs-data.inv_ext_hs_10;
    data.inv_ext_90=data.inv_ext-data.inv_ext_10;
    
    fig_inv_ext_10_hs=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_ext_hs_10(iiy:ify)./data.inv_ext_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs_90(iiy:ify)./data.inv_ext_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill External Inventors');
    xlim([iyear,fyear])

    %-------------
    % Adjustment
    %-------------

    %Skill adjusted (without entrant inventors)
    data.inv_tot_adj=data.inv_tot_hs+data.inv_tot_ls;
    data.inv_int_adj=data.inv_int_hs+data.inv_int_ls;
    data.inv_ext_adj=data.inv_ext_hs+data.inv_ext_ls;
    data.inv_ent_adj=data.inv_ent_hs+data.inv_ent_ls;

    data.inv_10_adj=data.inv_hs_10+data.inv_ls_10;
    data.inv_90_adj=data.inv_hs_90+data.inv_ls_90;

    %High skill inventors out of high-skill + low-skill
    fig_inv_tot_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_tot_hs(iiy:ify)./data.inv_tot_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Share of High Skill Inventors (Adjusted)'); ylim([0.49,0.51])
    xlim([iyear,fyear])

    % Share of high skill in internal and external (adjusted)
    fig_inv_int_ext_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs(iiy:ify)./data.inv_int_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs(iiy:ify)./data.inv_ext_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors (Adjusted)');
    xlim([iyear,fyear])

    % Share of high skill in internal, external, entrants (adjusted)
    fig_inv_int_ext_ent_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs(iiy:ify)./data.inv_int_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs(iiy:ify)./data.inv_ext_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_ent_hs(iiy:ify)./data.inv_ent_adj(iiy:ify),'--'); p3.Color=par.opt.green;
    lgd=legend([p1,p2,p3],{'Internal','External','Entrant'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors (Adjusted)');
    xlim([iyear,fyear])

    % Share of high skill in top 10 and bottom 90 (adjusted)
    fig_inv_10_90_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_hs_10(iiy:ify)./data.inv_10_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_hs_90(iiy:ify)./data.inv_90_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors (Adjusted)');
    xlim([iyear,fyear])

    %Share of hs internal in 10 and 90 (adjusted)
    data.inv_int_10_adj=data.inv_int_hs_10+data.inv_int_ls_10;
    data.inv_int_hs_90=data.inv_int_hs-data.inv_int_hs_10;
    data.inv_int_ls_90=data.inv_int_ls-data.inv_int_ls_10;
    data.inv_int_90_adj=data.inv_int_hs_90+data.inv_int_ls_90;
    
    fig_inv_int_10_90_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs_10(iiy:ify)./data.inv_int_10_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_int_hs_90(iiy:ify)./data.inv_int_90_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Internal Inventors');
    xlim([iyear,fyear])
    
    %Share of hs external in 10 and 90 (adjusted)
    data.inv_ext_10_adj=data.inv_ext_hs_10+data.inv_ext_ls_10;
    data.inv_ext_hs_90=data.inv_ext_hs-data.inv_ext_hs_10;
    data.inv_ext_ls_90=data.inv_ext_ls-data.inv_ext_ls_10;
    data.inv_ext_90_adj=data.inv_ext_hs_90+data.inv_ext_ls_90;
    
    fig_inv_ext_10_90_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_ext_hs_10(iiy:ify)./data.inv_ext_10_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs_90(iiy:ify)./data.inv_ext_90_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill External Inventors');
    xlim([iyear,fyear])

    %Share of hs internal and external  in 10  (adjusted)
    fig_inv_int_ext_10_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs_10(iiy:ify)./data.inv_int_10_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs_10(iiy:ify)./data.inv_ext_10_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors in Top 10 \%'); 
    xlim([iyear,fyear]); 
    %ylim([0.25,0.75])

    %Share of hs internal and external  in bottom 90 (adjusted)
    fig_inv_int_ext_90_hs_adj=figure; hold all;
    p1=plot(data.year(iiy:ify),data.inv_int_hs_90(iiy:ify)./data.inv_int_90_adj(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ext_hs_90(iiy:ify)./data.inv_ext_90_adj(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Share of High Skill Inventors in Bottom 90');
    xlim([iyear,fyear]); 
    
    %---------------
    % Quality
    %---------------
    
    %Patent quality top 10, bottom 90
    fig_q_10_90=figure; hold all;
    p1=plot(data.year(iiy:ify),data.q_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.q_90(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10','Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Patent Quality $Log(1+fcit3)$');
    xlim([iyear,fyear])
    
    %Patent quality incumbents and entrants
    fig_q_inc_ent=figure; hold all;
    p1=plot(data.year(iiy:ify),data.q_inc(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.q_ent(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Incumbents','Entrants'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Patent Quality $Log(1+fcit3)$');
    xlim([iyear,fyear])
    
    %Patent quality internal and external
    fig_q_int_ext=figure; hold all;
    p1=plot(data.year(iiy:ify),data.q_int(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.q_ext(iiy:ify),'--'); p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Internal','External'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Patent Quality $Log(1+fcit3)$');
    xlim([iyear,fyear])
    
    %Patent quality internal and external in top 10 and bottom 90
    fig_q_int_ext_10_90=figure; hold all;
    p1=plot(data.year(iiy:ify),data.q_int_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.q_ext_10(iiy:ify),'--'); p2.Color=par.opt.light_blue;
    p3=plot(data.year(iiy:ify),data.q_int_90(iiy:ify)); p3.Color=par.opt.maroon;
    p4=plot(data.year(iiy:ify),data.q_ext_90(iiy:ify),'--'); p4.Color=par.opt.maroon;
    lgd=legend([p1,p2,p3,p4],{'Internal Top 10','External Top 10','Internal Bottom 90','External Bottom 90'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Patent Quality $Log(1+fcit3)$');
    xlim([iyear,fyear])

    %----------------------
    % Inventor transitions 
    %----------------------

    %Probability low skill --> high skill
    fig_pr_ls_hs=figure; hold all;
    p0=plot(data.year(iiy:ify),data.inv_ls_np_hs(iiy:ify)./data.inv_ls_np(iiy:ify)); p0.Color='k';
    p1=plot(data.year(iiy:ify),data.inv_ls_np_hs_10(iiy:ify)./data.inv_ls_np_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ls_np_hs_90(iiy:ify)./data.inv_ls_np_90(iiy:ify)); p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_ls_np_hs_ent(iiy:ify)./data.inv_ls_np_ent(iiy:ify)); p3.Color=par.opt.green;
    lgd=legend([p0,p1,p2,p3],{'All','Top 10', 'Bottom 90', 'Entrants'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Probability low skill $\to$ high skill');
    xlim([iyear,fyear]); 

    %Probability high skill --> low skill
    fig_pr_hs_ls=figure; hold all;
    p0=plot(data.year(iiy:ify),data.inv_hs_np_ls(iiy:ify)./data.inv_hs_np(iiy:ify)); p0.Color='k';
    p1=plot(data.year(iiy:ify),data.inv_hs_np_ls_10(iiy:ify)./data.inv_hs_np_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_hs_np_ls_90(iiy:ify)./data.inv_hs_np_90(iiy:ify)); p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_hs_np_ls_ent(iiy:ify)./data.inv_hs_np_ent(iiy:ify)); p3.Color=par.opt.green;
    lgd=legend([p0,p1,p2,p3],{'All','Top 10', 'Bottom 90', 'Entrants'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Probability high skill $\to$ low skill');
    xlim([iyear,fyear]); 

    %Probability high skill --> out 
    fig_pr_hs_out=figure; hold all;
    p0=plot(data.year(iiy:ify),data.inv_hs_lp(iiy:ify)./data.inv_tot_hs(iiy:ify)); p0.Color='k';
    p1=plot(data.year(iiy:ify),data.inv_hs_lp_10(iiy:ify)./data.inv_hs_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_hs_lp_90(iiy:ify)./data.inv_hs_90(iiy:ify));  p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_hs_lp_ent(iiy:ify)./data.inv_ent_hs(iiy:ify));  p3.Color=par.opt.green;
    lgd=legend([p0,p1,p2,p3],{'All','Top 10', 'Bottom 90', 'Entrants'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Death Rate High Skill');
    xlim([iyear,fyear-dy]); 

    %Probability low skill --> out 
    fig_pr_ls_out=figure; hold all;
    p0=plot(data.year(iiy:ify),data.inv_ls_lp(iiy:ify)./data.inv_tot_ls(iiy:ify)); p0.Color='k';
    p1=plot(data.year(iiy:ify),data.inv_ls_lp_10(iiy:ify)./data.inv_ls_10(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.inv_ls_lp_90(iiy:ify)./data.inv_ls_90(iiy:ify));  p2.Color=par.opt.maroon;
    p3=plot(data.year(iiy:ify),data.inv_ls_lp_ent(iiy:ify)./data.inv_ent_ls(iiy:ify));  p3.Color=par.opt.green;
    lgd=legend([p0,p1,p2,p3],{'All','Top 10', 'Bottom 90', 'Entrants'}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Death Rate Low Skill');
    xlim([iyear,fyear-dy]); 

    %----------------------
    % Firm transitions 
    %----------------------

    %Probability of transition
    fig_pr_f_trans=figure; hold all;
    p1=plot(data.year(iiy:ify),data.firm_10_np_90(iiy:ify)./data.firm_10_np(iiy:ify)); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.firm_90_np_10(iiy:ify)./data.firm_90_np(iiy:ify));  p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10 $\to$ Bottom 90', 'Bottom 90 $\to$ Top 10 '}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Transition Probability');
    xlim([iyear,fyear-dy]); 

    %Probability of death
    fig_pr_f_out=figure; hold all;
    p1=plot(data.year(iiy:ify),data.firm_10_lp(iiy:ify)./(data.firm_10_np(iiy:ify)+data.firm_10_lp(iiy:ify))); p1.Color=par.opt.light_blue;
    p2=plot(data.year(iiy:ify),data.firm_90_lp(iiy:ify)./(data.firm_90_np(iiy:ify)+data.firm_90_lp(iiy:ify)));  p2.Color=par.opt.maroon;
    lgd=legend([p1,p2],{'Top 10', 'Bottom 90 '}); lgd.Location='best';  lgd.Color='w'; lgd.EdgeColor='none';
    xlabel('Year'); ylabel('Death Rate');
    xlim([iyear,fyear-dy]);

    %----------------------
    % Ratios
    %----------------------

    % Patents per inventor
    fig_pat_inv=figure; hold all;
    p1=plot(data.year(iiy:ify),data.pat_tot(iiy:ify)./data.inv_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Patents per Inventor');
    xlim([iyear,fyear])

    % Quality per inventor
    fig_lf_cit3_inv=figure; hold all;
    p1=plot(data.year(iiy:ify),data.logf_cit3(iiy:ify)./data.inv_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Quality per Inventor');
    xlim([iyear,fyear])

    % Value per inventor
    fig_lxi_inv=figure; hold all;
    p1=plot(data.year(iiy:ify),data.logxi(iiy:ify)./data.inv_tot(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Patent Value per Inventor');
    xlim([iyear,fyear])

    % Value per Quality
    fig_lf_cit3_lxi=figure; hold all;
    p1=plot(data.year(iiy:ify),data.logxi(iiy:ify)./data.logf_cit3(iiy:ify)); p1.Color=par.opt.light_blue;
    xlabel('Year'); ylabel('Patent Value per Citation');
    xlim([iyear,fyear])

    if save_fig==1
    %For Matlab 2020 onwards
        exportgraphics(fig_con,[par.opt.dir_fig  'data_mom_con' save_str '.pdf'])
        exportgraphics(fig_pat_inc,[par.opt.dir_fig  'data_mom_pat_inc' save_str '.pdf'])
        exportgraphics(fig_pat_int_ext,[par.opt.dir_fig  'data_mom_pat_int_ext' save_str '.pdf'])
        exportgraphics(fig_pat_int_10_90,[par.opt.dir_fig  'data_mom_pat_int_10_90' save_str '.pdf'])
        exportgraphics(fig_pat_inc_10,[par.opt.dir_fig  'data_mom_pat_inc_10' save_str '.pdf'])
        exportgraphics(fig_con_inv,[par.opt.dir_fig  'data_mom_con_inv' save_str '.pdf'])
        exportgraphics(fig_inv_inc,[par.opt.dir_fig  'data_mom_inv_inc' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext,[par.opt.dir_fig  'data_mom_inv_int_ext' save_str '.pdf'])
        exportgraphics(fig_inv_10_90_hs,[par.opt.dir_fig  'data_mom_inv_10_90_hs' save_str '.pdf'])
        exportgraphics(fig_inv_int_10_90,[par.opt.dir_fig  'data_mom_inv_int_10_90' save_str '.pdf'])
        exportgraphics(fig_inv_ext_10_90,[par.opt.dir_fig  'data_mom_inv_ext_10_90' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_10,[par.opt.dir_fig  'data_mom_inv_int_ext_10' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_90,[par.opt.dir_fig  'data_mom_int_ext_90' save_str '.pdf'])
        exportgraphics(fig_inv_tot_hs_ls_einv,[par.opt.dir_fig  'data_mom_inv_tot_hs_ls_einv' save_str '.pdf'])
        exportgraphics(fig_inv_hs_10,[par.opt.dir_fig  'data_mom_inv_hs_10' save_str '.pdf'])
        exportgraphics(fig_inv_inc_hs,[par.opt.dir_fig  'data_mom_inv_inc_hs' save_str '.pdf'])
        exportgraphics(fig_inv_hs_inc,[par.opt.dir_fig  'data_mom_inv_hs_inc' save_str '.pdf'])
        exportgraphics(fig_inv_ent_hs,[par.opt.dir_fig  'data_mom_inv_ent_hs' save_str '.pdf'])
        exportgraphics(fig_inv_hs_ent,[par.opt.dir_fig  'data_mom_inv_hs_ent' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_hs,[par.opt.dir_fig  'data_mom_inv_int_ext_hs' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_ent_hs,[par.opt.dir_fig  'data_mom_inv_int_ext_ent_hs' save_str '.pdf'])   
        exportgraphics(fig_inv_int_hs,[par.opt.dir_fig  'data_mom_inv_int_hs' save_str '.pdf'])
        exportgraphics(fig_inv_ext_hs,[par.opt.dir_fig  'data_mom_inv_ext_hs' save_str '.pdf'])
        exportgraphics(fig_inv_hs_int_ext_ent,[par.opt.dir_fig  'data_mom_inv_hs_int_ext_ent' save_str '.pdf'])
        exportgraphics(fig_inv_int_10_hs,[par.opt.dir_fig  'data_mom_inv_int_10_hs' save_str '.pdf'])
        exportgraphics(fig_inv_ext_10_hs,[par.opt.dir_fig  'data_mom_inv_ext_10_hs' save_str '.pdf'])
        exportgraphics(fig_inv_tot_hs_adj,[par.opt.dir_fig  'data_mom_inv_tot_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_hs_adj,[par.opt.dir_fig  'data_mom_inv_int_ext_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_ent_hs_adj,[par.opt.dir_fig  'data_mom_inv_int_ext_ent_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_10_90_hs_adj,[par.opt.dir_fig  'data_mom_inv_10_90_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_int_10_90_hs_adj,[par.opt.dir_fig  'data_mom_inv_int_10_90_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_ext_10_90_hs_adj,[par.opt.dir_fig  'data_mom_inv_ext_10_90_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_10_hs_adj,[par.opt.dir_fig  'data_mom_inv_int_ext_10_hs_adj' save_str '.pdf'])
        exportgraphics(fig_inv_int_ext_90_hs_adj,[par.opt.dir_fig  'data_mom_inv_int_ext_90_hs_adj' save_str '.pdf'])
        exportgraphics(fig_q_10_90,[par.opt.dir_fig  'data_mom_q_10_90' save_str '.pdf'])
        exportgraphics(fig_q_inc_ent,[par.opt.dir_fig  'data_mom_q_inc_ent' save_str '.pdf'])
        exportgraphics(fig_q_int_ext,[par.opt.dir_fig  'data_mom_q_int_ext' save_str '.pdf'])
        exportgraphics(fig_q_int_ext_10_90,[par.opt.dir_fig  'data_mom_q_int_ext_10_90' save_str '.pdf'])
        exportgraphics(fig_pr_ls_hs,[par.opt.dir_fig  'data_mom_pr_ls_hs' save_str '.pdf'])
        exportgraphics(fig_pr_hs_ls,[par.opt.dir_fig  'data_mom_pr_hs_ls' save_str '.pdf'])
        exportgraphics(fig_pr_hs_out,[par.opt.dir_fig  'data_mom_pr_hs_out' save_str '.pdf'])
        exportgraphics(fig_pr_ls_out,[par.opt.dir_fig  'data_mom_pr_ls_out' save_str '.pdf'])

        exportgraphics(fig_pr_f_trans,[par.opt.dir_fig  'data_mom_pr_f_trans' save_str '.pdf'])
        exportgraphics(fig_pr_f_out,[par.opt.dir_fig  'data_mom_pr_f_out' save_str '.pdf'])

        exportgraphics(fig_pat_inv,[par.opt.dir_fig  'data_mom_pat_inv' save_str '.pdf'])
        exportgraphics(fig_lf_cit3_inv,[par.opt.dir_fig  'data_mom_lf_cit3_inv' save_str '.pdf'])
        exportgraphics(fig_lxi_inv,[par.opt.dir_fig  'data_mom_lxi_inv' save_str '.pdf'])
        exportgraphics(fig_lf_cit3_lxi,[par.opt.dir_fig  'data_mom_lf_cit3_lxi' save_str '.pdf'])
       

    end
    
    
end

