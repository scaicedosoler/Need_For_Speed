%Make table of comparing fit of the two estimations of the model

function []=tab_fit_cmp_intq_B(dmom0,smm0,dmom1,smm1)

mom0=smm0.mom;
mom1=smm1.mom;

%Save table
save_tab=0;

%% Table of all parameters

%--------------------------------
% Parameters
%--------------------------------
  
%Unpack parameters
par0=smm0.par;
par1=smm1.par;

input.tableCaption='Parameters';
input.tableLabel='par';


%Resize option
resize=0;

if resize==1
resize_in='\resizebox{\textwidth}{!}{';
resize_out='}';
else
resize_in='';
resize_out='';
end

%Array stretch size
array_stretch=1.2;

% make table header lines:

header = '\begin{tabular}{clcc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}'];...
        ['\renewcommand{\arraystretch}{' num2str(array_stretch) '} ']; ...
        resize_in; header ;'\hline \hline'};

%Content of table
latex=[latex;...
       '\textbf{Parameter} & \textbf{Description} & \textbf{Value ' num2str(smm0.iyear) '-' num2str(smm0.iyear+5) ' } & \textbf{Value ' num2str(smm1.iyear) '-' num2str(smm1.iyear+5) '} \\';'\hline' ; ...
       '\multicolumn{3}{l}{\textit{Estimated using SMM}} &  \\ [2pt]'; ...
       
       %' \multicolumn{3}{l}{\textit{~~ --Levels (Both Periods)--}} &  \\'; ...
       ['$\chi$ & ' 'Productivity Incumbents & ' num2str(round(par0.chi,3)) '& ' num2str(round(par1.chi,3)) ' \\' ];...
       ['$\lambda$ & ' 'Scale Parameter of Quality & ' num2str(round(par0.lambda,3)) '& ' num2str(round(par1.lambda,3)) ' \\' ];...
       ['$\nu$ & ' 'Knowledge Diffusion for Entry & ' num2str(round(1-par0.nu,3)) '& ' num2str(round(1-par1.nu,3)) ' \\' ];... %S: Changed for easier interpretation (May 2024)
       ['$\chi_e$ & ' 'Productivity Entrants & ' num2str(round(par0.chi_e,3)) '& ' num2str(round(par1.chi_e,3)) ' \\' ];...
       ['$\lambda_e$ & ' 'Quality Step-Size Entrants & ' num2str(round(par0.lambda_e,3)) '& ' num2str(round(par1.lambda_e,3)) ' \\' ];...
       ['$\alpha_e$ & ' 'Entrant Innovation Elasticity of Labor & ' num2str(round(par0.alpha_e,3)) '& ' num2str(round(par1.alpha_e,3)) ' \\' ];...
       %['$\gamma_e$ & ' 'Elasticity of Entrant Quality & ' num2str(round(par0.gamma_e,3)) '& ' num2str(round(par1.gamma_e,3)) ' \\' ];...
       ['$\gamma_q$ & ' 'Elasticity of Incumbent Quality & ' num2str(round(par0.gamma_q,3)) '& ' num2str(round(par1.gamma_q,3)) ' \\ [2pt]' ];...
       %' \multicolumn{3}{l}{\textit{~~ --Changes (Only 2010)--}} &  \\'; ...
       ['$L_{I}$ & ' 'Supply of Inventors & ' num2str(round(par0.L_I,3)) '& ' num2str(round(par1.L_I,3)) ' \\' ];...
       ['$\chi_b$ & ' 'Additional Benefit & ' num2str(round(par0.chi_b,2)) '& ' num2str(round(par1.chi_b,2)) ' \\' ];...  
       '\multicolumn{3}{l}{\textit{Matched From Regressions}} &  \\ [2pt]'; ...
       ['$\alpha_x$ & ' 'Speed Elasticity of Labor & ' num2str(round(par0.alpha_x,3)) '& ' num2str(round(par1.alpha_x,3)) ' \\' ];...
       ['$\alpha_q$ & ' 'Quality Elasticity of Labor & ' num2str(round(par0.alpha_q,3)) '& ' num2str(round(par1.alpha_q,3)) ' \\ [2pt]' ];... 
       %'\multicolumn{3}{l}{\textit{Estimated from Policy Shock}} &  \\'; ...
       '\multicolumn{3}{l}{\textit{Externally Calibrated and Normalized Parameters}} &  \\'; ...
       
       ['$\rho$ & ' 'Discount rate & ' num2str(round(par0.rho,2)) '& ' num2str(round(par1.rho,2)) ' \\' ];...
       ['$\beta$ & ' 'Elasticity of Labor of Production Workers & ' num2str(round(par0.beta,2)) '& ' num2str(round(par1.beta,2)) ' \\' ];...
       ['$L_p$ & ' 'Supply of Production Labor & ' num2str(round(par0.Lp,1)) '& ' num2str(round(par1.Lp,1)) ' \\' ];...
       % ['$\bar \gamma_q$ & ' 'Elasticity of Average Quality for Incumbents & ' num2str(round(par0.bargamma_q,2)) '& ' num2str(round(par1.bargamma_q,2)) ' \\' ];...
       % ['$\bar \gamma_e$ & ' 'Elasticity of Average Quality for Entrants & ' num2str(round(par0.bargamma_e,2)) '& ' num2str(round(par1.bargamma_e,2)) ' \\' ];...
       % '\multicolumn{3}{l}{\textit{Implied by Balanced Growth}} &  \\'; ...
       % ['$\gamma_x$ & ' 'Elasticity of Quality on Speed for Incumbents & ' num2str(round(par0.gamma_x,2)) '& ' num2str(round(par1.gamma_x,2)) ' \\' ];...
       % ['$\bar \gamma_x$ & ' 'Elasticity of Average Quality on Speed for Incumbents & ' num2str(round(par0.bargamma_x,2)) '& ' num2str(round(par1.bargamma_x,2)) ' \\' ];...
       % ['$\xi_e$ & ' 'Elasticity of Quality on Speed for Entrants & ' num2str(round(par0.xi_e,2)) '& ' num2str(round(par1.xi_e,2)) ' \\' ];...
       % ['$\bar \xi_e$ & ' 'Elasticity of Average Quality on Speed for Entrants & ' num2str(round(par0.barxi_e,2)) '& ' num2str(round(par1.barxi_e,2)) ' \\' ];...
       ];
   
footer = {'\hline \hline';
    ['\multicolumn{3}{l}{\scriptsize \textit{Notes}: Estimated parameters for the period '  num2str(smm0.iyear) '-' num2str(smm0.iyear+5) ' and '  num2str(smm1.iyear) '-' num2str(smm1.iyear+5) '.} ']
'\end{tabular}'; resize_out ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('%%---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));

fprintf('\n')
disp('\FloatBarrier')


%%% Save table

if par0.opt.save_tab==1
    %Store it in .tex document
    dlmcell([par.opt.dir_tab 'par_cmp.tex'],latex)
end

%% Moments

%Missing moments

if ~isfield(mom0,"r_lxi_lq")
    mom0.r_lxi_lq=nan;
    mom1.r_lxi_lq=nan;
end

%--------------------------------
% Moments Fit
%--------------------------------

input.tableCaption='Moments Fit';
input.tableLabel='mom_fit';
      
resize=1;

if resize==1
resize_in='\resizebox{\textwidth}{!}{';
resize_out='}';
else
resize_in='';
resize_out='';
end

%Array stretch size
array_stretch=1.2;

% make table header lines:

header = '\begin{tabular}{lcc|cc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}'];...
        ['\renewcommand{\arraystretch}{' num2str(array_stretch) '} ']; ...
        resize_in; header ;'\hline \hline'};

if ~isfield(dmom1, 'rb_pat_val_sales')

%For older versions of calibration
latex=[latex;...
       ['\textbf{Moment} & \textbf{Data ' num2str(smm0.iyear) '} & \textbf{Model ' num2str(smm0.iyear) '} & \textbf{Data ' num2str(smm1.iyear) '} & \textbf{Model ' num2str(smm1.iyear) '}  \\'];'\hline'; ...
       ' \multicolumn{2}{l}{\textit{--Levels (Both Periods)--}} & & & \\'; ...
       [' Average Growth Rate & '  num2str(100*round(dmom0.g,4)) '\% &' num2str(100*round(mom0.g,4)) '\% & '  num2str(100*round(dmom1.g,4)) '\% &' num2str(100*round(mom1.g,4)) '\%  \\' ];...    
       ['Inventors in Entrant Firms & '  num2str(100*round(dmom0.inv_ent,3)) '\% &' num2str(100*round(mom0.inv_ent,3)) '\% & '  num2str(100*round(dmom1.inv_ent,3)) '\% &' num2str(100*round(mom1.inv_ent,3)) '\% \\' ];...
       ['Share of Entrant Innovations & '  num2str(100*round(dmom0.pat_ent,3)) '\% &' num2str(100*round(mom0.pat_ent,3)) '\% & '  num2str(100*round(dmom1.pat_ent,3)) '\% &' num2str(100*round(mom1.pat_ent,3)) '\%  \\' ];...
       ['Patent Quality: $Q_{ent}/Q_{inc}$ & '  num2str(round(dmom0.q_ent_inc,3))  '&' num2str(round(mom0.q_ent_inc,3)) ' & '  num2str(round(dmom1.q_ent_inc,3))  ' &' num2str(round(mom1.q_ent_inc,3)) ' \\' ];...
       ['Inventors in Top 10 & '   num2str(100*round(dmom0.inv_10,3))  '\% &' num2str(100*round(mom0.inv_10,3)) '\% & '  num2str(100*round(dmom1.inv_10,3))  '\% &' num2str(100*round(mom1.inv_10,3)) '\%  \\' ];...
       ['Share of Innovations in Top 10& '   num2str(100*round(dmom0.pat_10,3))  '\% &' num2str(100*round(mom0.pat_10,3)) '\% & '  num2str(100*round(dmom1.pat_10,3))  '\% &' num2str(100*round(mom1.pat_10,3)) '\%  \\' ];...
       ['Patent Quality: Top 10 Relative to Bottom 90 & '  num2str(round(dmom0.q_10_90,3))  '&' num2str(round(mom0.q_10_90,3)) '& '  num2str(round(dmom1.q_10_90,3))  '&' num2str(round(mom1.q_10_90,3)) ' \\' ];...
       ' \multicolumn{2}{l}{\textit{--Changes (Only 2010)--}} & & &  \\'; ...
       ['Change in Speed: ($x/x_{base~year}$) & - &- & ' num2str(round(dmom1.pat_f_rb,2)) ' &' num2str(round(mom1.pat_f_rb,2)) '\\' ];...
       ['Change in Private Value/Quality: ($\log \left(V/Q \right )$) & - &- & '  num2str(round(dmom1.rb_xi_lq,2)) '&' num2str(round(mom1.rb_xi_lq,2)) '\\' ];...
       ['Change in Inventors per Patent  & - &- & '  num2str(round(dmom1.rb_inv_pat,2)) '&' num2str(round(mom1.rb_inv_pat,2)) '\\' ];...
        ];
else

%Most recent moments: May 2024
latex=[latex;...
       ['\textbf{Moment} & \textbf{Data ' num2str(smm0.iyear) '} & \textbf{Model ' num2str(smm0.iyear) '} & \textbf{Data ' num2str(smm1.iyear) '} & \textbf{Model ' num2str(smm1.iyear) '}  \\'];'\hline'; ...
       %' \multicolumn{2}{l}{\textit{--Levels (Both Periods)--}} & & & \\'; ...
       [' Average Growth Rate & '  num2str(100*round(dmom0.g,4)) '\% &' num2str(100*round(mom0.g,4)) '\% & '  num2str(100*round(dmom1.g,4)) '\% &' num2str(100*round(mom1.g,4)) '\%  \\' ];...    
       ['Inventors in Entrant Firms & '  num2str(100*round(dmom0.inv_ent,3)) '\% &' num2str(100*round(mom0.inv_ent,3)) '\% & '  num2str(100*round(dmom1.inv_ent,3)) '\% &' num2str(100*round(mom1.inv_ent,3)) '\% \\' ];...
       ['Share of Entrant Innovations & '  num2str(100*round(dmom0.pat_ent,3)) '\% &' num2str(100*round(mom0.pat_ent,3)) '\% & '  num2str(100*round(dmom1.pat_ent,3)) '\% &' num2str(100*round(mom1.pat_ent,3)) '\%  \\' ];...
       ['Patent Quality: $Q_{ent}/Q_{inc}$ & '  num2str(round(dmom0.q_ent_inc,3))  '&' num2str(round(mom0.q_ent_inc,3)) ' & '  num2str(round(dmom1.q_ent_inc,3))  ' &' num2str(round(mom1.q_ent_inc,3)) ' \\' ];...
       ['Inventors in Top 10 & '   num2str(100*round(dmom0.inv_10,3))  '\% &' num2str(100*round(mom0.inv_10,3)) '\% & '  num2str(100*round(dmom1.inv_10,3))  '\% &' num2str(100*round(mom1.inv_10,3)) '\%  \\' ];...
       ['Share of Innovations in Top 10& '   num2str(100*round(dmom0.pat_10,3))  '\% &' num2str(100*round(mom0.pat_10,3)) '\% & '  num2str(100*round(dmom1.pat_10,3))  '\% &' num2str(100*round(mom1.pat_10,3)) '\%  \\' ];...
       ['Patent Quality: Top 10 Relative to Bottom 90 & '  num2str(round(dmom0.q_10_90,3))  '&' num2str(round(mom0.q_10_90,3)) '& '  num2str(round(dmom1.q_10_90,3))  '&' num2str(round(mom1.q_10_90,3)) ' \\' ];...
       ['Patents per Unique Inventor & '  num2str(round(dmom0.x_l,3))  '&' num2str(round(mom0.x_l,3)) '& '  num2str(round(dmom1.x_l,3))  '&' num2str(round(mom1.x_l,3)) ' \\' ];...
        ' \multicolumn{2}{l}{\textit{--Changes (Relative to Baseline Year)--}} & & &  \\'; ...
       ['Change in Speed: ($x/x_{base~year}$) & - &- & ' num2str(round(dmom1.pat_f_rb,2)) ' &' num2str(round(mom1.pat_f_rb,2)) '\\' ];...
       ['Change in  Patents Value/Sales & - &- & '  num2str(round(dmom1.rb_pat_val_sales,2))  '&' num2str(round(mom1.rb_pat_val_sales,2)) ' \\' ];...'
      
       % ['Change in Private Value/Quality: ($\log \left(V/Q \right )$) & - &- & '  num2str(round(dmom1.rb_xi_lq,2)) '&' num2str(round(mom1.rb_xi_lq,2)) '\\' ];...
       % ['Change in Inventors per Patent  & - &- & '  num2str(round(dmom1.rb_inv_pat,2)) '&' num2str(round(mom1.rb_inv_pat,2)) '\\' ];...
        ];
end

       
footer = {'\hline \hline';
    ['\multicolumn{3}{l}{\scriptsize \textit{Notes}: Data targeted moments and model fit for the period '  num2str(smm0.iyear) '-' num2str(smm0.iyear+5) ' and '  num2str(smm1.iyear) '-' num2str(smm1.iyear+5) '.} ']
'\end{tabular}'; resize_out ;'\end{table}'};

latex = [latex;footer];

fprintf('\n')
fprintf('%%---------------------------------------------------------------------')
fprintf('\n')

%Print table to copy
disp(char(latex));

fprintf('\n')
disp('\FloatBarrier')


%%% Save table

if save_tab==1
    %Store it in .tex document
    dlmcell([par.opt.dir_tab 'fit_cmp.tex'],latex)
end
