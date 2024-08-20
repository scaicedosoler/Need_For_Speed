%Growth decomposition

function g_deco_intq_B(eq0,eq1,par0,par1)

%--------------------------------
% Growth decomposition
%--------------------------------
eq0.Qe=eq0.dqe;
eq1.Qe=eq1.dqe;

eq0.Q=eq0.dq;
eq1.Q=eq1.dq;

%Growth entrants and Incumbents
eq0.g_ent=eq0.xe*eq0.Qe; eq0.g_inc=eq0.x*eq0.Q;
eq1.g_ent=eq1.xe*eq1.Qe; eq1.g_inc=eq1.x*eq1.Q;

%Change in growth
eq1.Dg_ent=(eq1.xe*eq1.Qe)/(eq0.xe*eq0.Qe)-1;
eq1.Dg_inc=(eq1.x*eq1.Q)/(eq0.x*eq0.Q)-1;
eq1.Dg=eq1.g/eq0.g-1;

%Change in speed and quality
eq1.Dx_ent=eq1.xe/eq0.xe-1;
eq1.Dx_inc=eq1.x/eq0.x-1;
eq1.Dq_ent=eq1.Qe/eq0.Qe-1;
eq1.Dq_inc=eq1.Q/eq0.Q-1;

eq1.Dxq_ent=eq1.Dx_ent*eq1.Dq_ent;
eq1.Dxq_inc=eq1.Dx_inc*eq1.Dq_inc;


eq1.Dx=(eq0.xe/(eq0.xe+eq0.x))*eq1.Dx_ent+(eq0.x/(eq0.xe+eq0.x))*eq1.Dx_inc;
eq1.Dq=(eq0.Qe/(eq0.Qe+eq0.Q))*eq1.Dq_ent+(eq0.Q/(eq0.Qe+eq0.Q))*eq1.Dq_inc;

eq1.Dxq=eq1.Dg-eq1.Dx-eq1.Dq;%eq1.Dx*eq1.Dq;

%Table

input.tableCaption='Growth Decomposition';
input.tableLabel='g_deco';

resize=0;

if resize==1
    resize_in='\resizebox{\textwidth}{!}{';
    resize_out='}';
else
    resize_in='';
    resize_out='';
end

% make table header lines:

header = '\begin{tabular}{lcccccc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; ...
    header ;'\hline \hline'};

% latex = {'\begin{table}[h!]';'\centering'; ...
%     resize_in; header ;'\hline \hline'};

%latex = {'\begin{table}[h!]';'\centering'; header ;'\hline \hline'};

%Content of table
latex=[latex;...
    '& \multicolumn{2}{c}{\textbf{Growth Rate}} & \multicolumn{4}{c}{\textbf{Percentage Change}} \\';...
    '\textbf{Firms} & \textit{' num2str(par0.iyear) '}  & \textit{' num2str(par1.iyear) '}   & \textit{Speed} & \textit{Quality} & \textit{Speed $\times$ Quality}& \textit{Total}  \\';'\hline'; ...
    ['Entrant & '  num2str(100*round(eq0.g_ent,4)) '\% & ' num2str(100*round(eq1.g_ent,4))   '\% & ' num2str(100*round(eq1.Dx_ent,3)) '\% & ' num2str(100*round(eq1.Dq_ent,3)) '\% & ' num2str(100*round(eq1.Dxq_ent,3)) '\% & ' num2str(100*round(eq1.Dg_ent,3)) '\% \\' ];...
    ['Incumbent & '  num2str(100*round(eq0.g_inc,4)) '\% & ' num2str(100*round(eq1.g_inc,4))  '\% & ' num2str(100*round(eq1.Dx_inc,3)) '\% & ' num2str(100*round(eq1.Dq_inc,3)) '\% & ' num2str(100*round(eq1.Dxq_inc,3)) '\% & ' num2str(100*round(eq1.Dg_inc,3))  '\% \\' ];...
    ['Total & '  num2str(100*round(eq0.g,4)) '\% & ' num2str(100*round(eq1.g,4))  '\% & ' num2str(100*round(eq1.Dx,3)) '\% & ' num2str(100*round(eq1.Dq,3)) '\% & ' num2str(100*round(eq1.Dxq,3))  '\% & ' num2str(100*round(eq1.Dg,3))  '\% \\' ];...
    ];

%footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};
footer = {'\hline \hline';'\end{tabular}' ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('%%---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));

fprintf('\n')
disp('\FloatBarrier')

%--------------------------------
% Labor allocation
%--------------------------------

%Share of labor
eq0.sle=eq0.le/par0.L_I; eq1.sle=eq1.le/par1.L_I;
eq0.slx=eq0.lx/par0.L_I; eq1.slx=eq1.lx/par1.L_I;
eq0.slq=eq0.lq/par0.L_I; eq1.slq=eq1.lq/par1.L_I;

%Table

input.tableCaption='Allocation of Inventors';
input.tableLabel='inv_allocation';

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
    header ;'\hline \hline'};

% latex = {'\begin{table}[h!]';'\centering'; ...
%     resize_in; header ;'\hline \hline'};

%latex = {'\begin{table}[h!]';'\centering'; header ;'\hline \hline'};

%Content of table
latex=[latex;...
    '\textbf{Inventors} & \textbf{' num2str(par0.iyear) '}  & \textbf{' num2str(par1.iyear) '}   \\' ; ...
    ['Speed (Incumbents) & '  num2str(100*round(eq0.slx,3)) '\% & ' num2str(100*round(eq1.slx,3))   '\%  \\' ];...
    ['Quality (Incumbents) & '  num2str(100*round(eq0.slq,3)) '\% & ' num2str(100*round(eq1.slq,3))   '\%  \\' ];...
    ['Entrant Firms & '  num2str(100*round(eq0.sle,3)) '\% & ' num2str(100*round(eq1.sle,3))   '\%  \\' ];...
    ];

%footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};
footer = {'\hline \hline';'\end{tabular}' ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('%%---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));

fprintf('\n')
disp('\FloatBarrier')

end