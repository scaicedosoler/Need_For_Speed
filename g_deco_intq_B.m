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

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; '\renewcommand{\arraystretch}{1.2}' ; ...
    header ;'\hline \hline'};

% latex = {'\begin{table}[h!]';'\centering'; ...
%     resize_in; header ;'\hline \hline'};

%latex = {'\begin{table}[h!]';'\centering'; header ;'\hline \hline'};

%Content of table
latex=[latex;...
    '& \multicolumn{2}{c}{\textbf{Growth Rate}} & \multicolumn{4}{c}{\textbf{Percentage Change}} \\';...
    '\textbf{Firms} & \textit{' num2str(par0.iyear) '}  & \textit{' num2str(par1.iyear) '}   & \textit{Speed} & \textit{Quality} & \textit{Speed $\times$ Quality}& \textit{Total}  \\';'\hline'; ...
    ['Incumbent & '  num2str(100*round(eq0.g_inc,4)) '\% & ' num2str(100*round(eq1.g_inc,4))  '\% & ' num2str(100*round(eq1.Dx_inc,3)) '\% & ' num2str(100*round(eq1.Dq_inc,3)) '\% & ' num2str(100*round(eq1.Dxq_inc,3)) '\% & ' num2str(100*round(eq1.Dg_inc,3))  '\% \\' ];...
    ['Entrant & '  num2str(100*round(eq0.g_ent,4)) '\% & ' num2str(100*round(eq1.g_ent,4))   '\% & ' num2str(100*round(eq1.Dx_ent,3)) '\% & ' num2str(100*round(eq1.Dq_ent,3)) '\% & ' num2str(100*round(eq1.Dxq_ent,3)) '\% & ' num2str(100*round(eq1.Dg_ent,3)) '\% \\' ];...
    ['Total & '  num2str(100*round(eq0.g,4)) '\% & ' num2str(100*round(eq1.g,4))  '\% & ' num2str(100*round(eq1.Dx,3)) '\% & ' num2str(100*round(eq1.Dq,3)) '\% & ' num2str(100*round(eq1.Dxq,3))  '\% & ' num2str(100*round(eq1.Dg,3))  '\% \\' ];...
    ];

%footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};
footer = {'\hline \hline';...
    '\multicolumn{7}{l}{\scriptsize \textit{Notes}: Growth rate decomposition into percentage changes for incumbents and entrants.}  '; ...
    '\end{tabular}' ;'\end{table}'};

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

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; '\renewcommand{\arraystretch}{1.2}' ; ...
    header ;'\hline \hline'};

% latex = {'\begin{table}[h!]';'\centering'; ...
%     resize_in; header ;'\hline \hline'};

%latex = {'\begin{table}[h!]';'\centering'; header ;'\hline \hline'};

%Content of table
latex=[latex;...
    '\textbf{Inventors} & \textbf{' num2str(par0.iyear) '}  & \textbf{' num2str(par1.iyear) '}   \\' ; '\hline'; ...
    ['Speed (Incumbents) & '  num2str(100*round(eq0.slx,3)) '\% & ' num2str(100*round(eq1.slx,3))   '\%  \\' ];...
    ['Quality (Incumbents) & '  num2str(100*round(eq0.slq,3)) '\% & ' num2str(100*round(eq1.slq,3))   '\%  \\' ];...
    ['Entrant Firms & '  num2str(100*round(eq0.sle,3)) '\% & ' num2str(100*round(eq1.sle,3))   '\%  \\' ];...
    ];

%footer = {'\hline \hline';'\end{tabular}'; resize_out ;'\end{table}'};
footer = {'\hline \hline';...
    '\multicolumn{2}{l}{\scriptsize \textit{Notes}: Allocation of inventors across speed and quality in 1980-1985 and 2010-2015.}   '; ...
    '\end{tabular}' ;'\end{table}'};

latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('%%---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));

fprintf('\n')
disp('\FloatBarrier')

%--------------------------------------------------------------
% Growth Decomposition: change in parameters and variables
%--------------------------------------------------------------
linc1=eq1.lx+eq1.lq;
linc0=eq0.lx+eq0.lq;

%Log decomposition--not so accurate approximation
Dlg=log(eq1.g)-log(eq0.g);

%Entry
we0=eq0.xe*eq0.dqe/eq0.g;
ge1=eq1.xe*eq1.dqe;
ge0=eq0.xe*eq0.dqe;
Dlge=log(ge1)-log(ge0);

%Incumbents
winc0=eq0.x*eq0.dq/eq0.g;
ginc1=eq1.x*eq1.dq;
ginc0=eq0.x*eq0.dq;
Dlginc=log(ginc1)-log(ginc0);


%Quality
DlQ=log(eq1.dq)-log(eq0.dq);

Dllambda=log(par1.lambda)-log(par0.lambda);
Dllq=log(eq1.lq)-log(eq0.lq);
Dalpha_q=par1.alpha_q-par0.alpha_q;

Dlq_d1=Dllambda;
Dlq_d2=par1.alpha_q*Dllq;
Dlq_d3=Dalpha_q*log(eq0.lq);

Dlinc=log(linc1)-log(linc0);
Dllq_linc=Dllq-Dlinc;


%Speed
Dlx=log(eq1.x)-log(eq0.x);

Dlchi=log(par1.chi)-log(par0.chi);
Dllx=log(eq1.lx)-log(eq0.lx);
Dalpha_x=par1.alpha_x-par0.alpha_x;

Dlx_d1=Dlchi;
Dlx_d2=par1.alpha_x*Dllx;
Dlx_d3=Dalpha_x*log(eq0.lx);

Dlinc=log(linc1)-log(linc0);
Dllx_linc=Dllx-Dlinc;


%Table

input.tableCaption='Speed and Quality Decomposition';
input.tableLabel='xq_deco';

resize=1;

if resize==1
    resize_in='\resizebox{\textwidth}{!}{';
    resize_out='}';
else
    resize_in='';
    resize_out='';
end

% make table header lines:

%Array stretch size
array_stretch=1.2;

% make table header lines:

header = '\begin{tabular}{lcccccc}';

% latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; ...
%     header ;'\hline \hline'};

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}'];...
        ['\renewcommand{\arraystretch}{' num2str(array_stretch) '} ']; ...
        resize_in; header ;'\hline \hline'};

%latex = {'\begin{table}[h!]';'\centering'; header ;'\hline \hline'};

%Content of table
latex=[latex;...
    ' & Productivity   & Labor Quantity  & Labor Allocation & Elasticity   & Total  \\ [-3pt]' ; ...
    '& \rowfont{\scriptsize} $\Delta \log(\chi) \vee \Delta \log(\lambda)  $ \rowfont{\scriptsize} & \rowfont{\scriptsize} $\alpha_i \Delta \log(l_{inc}) $  & \rowfont{\scriptsize} $\alpha_i \Delta \log(l_i/l_{inc}) $ & \rowfont{\scriptsize} $\Delta \alpha_i \log(l_i)$   & \\   ' ; ...
    '\hline' ;
    [' \textit{Speed} & '  num2str(round(Dlchi,3)) ' & ' num2str(round(par1.alpha_x*Dlinc,3)) ' & ' num2str(round(par1.alpha_x*Dllx_linc,3)) '& ' num2str(round(Dalpha_x*log(eq0.lx),3)) '& ' num2str(round(Dlx,3))   ' \\' ];...
    [' \textit{Quality}& '  num2str(round(Dllambda,3)) ' & ' num2str(round(par1.alpha_q*Dlinc,3)) ' & ' num2str(round(par1.alpha_q*Dllq_linc,3)) '& ' num2str(round(Dalpha_q*log(eq0.lq),3)) '& ' num2str(round(DlQ,3))   ' \\' ];...

    ];

footer = {'\hline \hline'; '\multicolumn{5}{l}{\scriptsize \textit{Notes}: Speed and quality decomposition into changes in productivity, labor and elasticities.} ';'\end{tabular}'; resize_out ;'\end{table}'};
%footer = {'\hline \hline';'\end{tabular}' ;'\end{table}'};



latex = [latex;footer];

%Print table to copy
fprintf('\n')
fprintf('%%---------------------------------------------------------------------')
fprintf('\n')

disp(char(latex));

fprintf('\n')
disp('\FloatBarrier')

end