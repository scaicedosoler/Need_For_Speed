%Make table of parameters and variables

function []=tab_fun(par,eq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input.tableCaption='Parameters';
input.tableLabel='par';

resize=1;

if resize==1
resize_in='\resizebox{\textwidth}{!}{';
resize_out='}';
else
resize_in='';
resize_out='';
end

% make table header lines:
header = '\begin{tabular}{clc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; ...
    resize_in; header ;'\hline \hline'};

%Content of table
latex=[latex;...
       '\textbf{Parameter} & \textbf{Description} & \textbf{Value} \\';'\hline' ; ...
       '\multicolumn{1}{l}{\textbf{Labor Supply}} & & \\'; ...
       ['$L_{I}$ & ' 'Supply of Inventors & ' num2str(round(par.L_I,2)) ' \\' ];...

        '\multicolumn{1}{l}{\textbf{Arrival Rates}} & & \\'; ...
       ['$\alpha_x$ & ' 'Incumbent Innovation Elasticity of Labor & ' num2str(round(par.alpha_x,4)) ' \\' ];...
       ['$\chi$ & ' 'Radical Innovation Productivity for Low Type Firms & ' num2str(round(par.chi,2)) ' \\' ];...
       ['$\chi_e$ & ' 'Entrants Innovation Productivity  & ' num2str(round(par.chi_e,2)) ' \\' ];...
       '\multicolumn{1}{l}{\textbf{Quality}} & & \\'; ...
       ['$\lambda$ & ' 'Radical Innovation Step Size for High Type Firms & ' num2str(round(par.lambda,3)) ' \\' ];...
       ['$\alpha_q$ & ' 'Incumbent Quality Elasticity of Labor & ' num2str(round(par.alpha_q,4)) ' \\' ];...
       ['$\lambda_e$ & ' 'Quality Step-Size for Entrants & ' num2str(round(par.lambda_e,3)) ' \\' ] ...
     
       ];
       
footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('\n ---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));


%%% Save table

if par.opt.save_tab==1
    %Store it in .tex document
    dlmcell([par.opt.dir_tab 'par.tex'],latex)
end

fprintf('\n')
disp('\FloatBarrier')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table of equilibrium variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Labor
eq.slx=eq.Phi*eq.lx/par.L_I;
eq.slq=eq.Phi*eq.lq/par.L_I;
eq.sle=eq.le/par.L_I;
eq.lqx=eq.lq/(eq.lq+eq.lx);

eq.w=eq.wq*par.qbar;

%Innovation
eq.sx=eq.Phi*eq.x/(eq.Phi*eq.x+eq.xe);

input.tableCaption='Equilibrium Variables';
input.tableLabel='eq_var';
      
resize=0;

if resize==1
resize_in='\resizebox{\textwidth}{!}{';
resize_out='}';
else
resize_in='';
resize_out='';
end


% make table header lines:
header = '\begin{tabular}{clc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; ...
    resize_in; header ;'\hline \hline'};

%Content of table
latex=[latex;...
       '\textbf{Variable} & \textbf{Description} & \textbf{Value} \\';'\hline' ; ...
       ['$g$ & ' 'Growth Rate & ' num2str(100*round(eq.g,4)) '\% \\' ];...
       '\multicolumn{1}{l}{\textbf{Labor Variables}} & \\'; ...
       ['$s_{lx}$ & ' 'Share of Labor in Incumbent Arrival Rate & ' num2str(100*round(eq.slx,2)) '\% \\' ];...
       ['$s_{lq}$ & ' 'Share of Labor in Incumbent Quality & ' num2str(100*round(eq.slq,2)) '\% \\' ];...
       ['$s_{le}$ & ' 'Share of Labor in Entrants & ' num2str(100*round(eq.sle,2)) '\% \\' ];...
       ['$\frac{l_q}{l_x+l_q}$ & ' 'Share of Labor for in Quality & ' num2str(100*round(eq.lqx,2)) '\% \\' ];...
       
       '\multicolumn{1}{l}{\textbf{Innovation and Quality}} & \\'; ...
       ['$s_x$ & ' 'Share of Incumbents in Innovation & ' num2str(100*round(eq.sx,2)) ' \% \\' ];...
       ['$\Phi $ & ' 'Measure of Firms & ' num2str(round(eq.Phi,3)) ' \\' ];...
       ['$\Phi x$ & ' 'Incumbents Innovation Rate & ' num2str(round(eq.Phi*eq.x,4)) ' \\' ];...
       ['$x_e$ & ' 'Entrant Innovation Rate & ' num2str(round(eq.xe,4)) ' \\' ];...
       ['$Q$ & ' 'Quality of Innovation (Step-Size) Incumbents & ' num2str(round(eq.Q,3)) ' \\' ];...
       '\multicolumn{1}{l}{\textbf{Wages}} & \\'; ...
       ['$w$ & ' 'Wages& ' num2str(round(eq.w,3)) ' \\' ]...
       ];
       
footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('\n ---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));


%%% Save table

if par.opt.save_tab==1
    %Store it in .tex document
    dlmcell([par.opt.dir_tab 'var.tex'],latex)
end

fprintf('\n')
disp('\FloatBarrier')

%% Table of moments

%{
tab_mom=1;

if tab_mom==1

input.tableCaption='Moments to Match';
input.tableLabel='mom_match';


resize=0;

if resize==1
resize_in='\resizebox{\textwidth}{!}{';
resize_out='}';
else
resize_in='';
resize_out='';
end

% make table header lines:

header = '\begin{tabular}{lcc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; ...
    resize_in; header ;'\hline \hline'};

%Content of table
latex=[latex;...
        ['\textbf{Moment} & \textbf{Model}  \\'];'\hline'; ...
        ['1. Average growth rate &' num2str(round(eq.g,3)) ' \\' ];...
        '\textbf{Share of Inventors}  & \\ ';...
        ['2. Inventors in Top 10 &' num2str(round(eq.sl_top,3)) ' \\' ];...
        ['3. Inventors in Entrant Firms &' num2str(round(eq.sle,3)) ' \\' ];...
        ['4. High Skill Inventors in Top 10 &' num2str(round(eq.slh_top,3)) ' \\' ];...
        ['5. High Skill Inventors in Entrant Firms &' num2str(round(eq.sleh,3)) ' \\' ];...
        '\textbf{Share of Innovations}  &  \\ ';...
        ['6. Innovations by Top 10 &' num2str(round(eq.s_top,3)) ' \\' ];...
        ['7. Share of Entrant Innovations &' num2str(round(eq.se,3)) ' \\' ];...
        '\textbf{Patent Quality}  &  \\ ';...
        ['8. Patent Quality Top 10 Relative to Average &' num2str(round(eq.qr_top_qbar_adj,3)) ' \\' ];...
        ['9. Patent Quality in Bottom 90  Relative to Average &' num2str(round(eq.qr_bot_qbar_adj,3)) ' \\' ];...
        ];

footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('\n ---------------------------------------------------------------------')
fprintf('\n')

fprintf(char(latex));


%%% Save table

if par.opt.save_tab==1
    %Store it in .tex document
    dlmcell([par.opt.dir_tab 'mom_tab.tex'],latex)
end

fprintf('\n')
fprintf('\FloatBarrier')

%}

end


