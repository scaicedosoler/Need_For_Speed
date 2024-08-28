%Counterfactual: Changing alpha_q and alpha_x
function cf_alpha_q_alpha_x(par,eq,par0,eq0)

if ~isfield(par,'byear')
    par.byear=2010;
    par.cyear=1980;
end

%-------------------------------------
% Reallocating labor
%-------------------------------------

%1. Changing quantity of labor
eqc1.lx=eq0.lx;
eqc1.lq=eq0.lq;
eqc1.le=eq0.le;

% Arrival rate
eqc1.x=par.chi*eqc1.lx^par.alpha_x;

%Quality
eqc1.dq=par.lambda*eqc1.lq^par.alpha_q;

%Entry
eqc1.xe=par.chi_e*eqc1.le^par.alpha_e;
eqc1.dqe=par.lambda_e;

%Growth
eqc1.g=eqc1.xe*eqc1.dqe+eqc1.x*eqc1.dq;




%2. Changing proportion of speed and quality of labor

eqc2.p_lx=eq0.lx/(eq0.lx+eq0.lq);

eqc2.lx=(eq.lx+eq.lq)*eqc2.p_lx;
eqc2.lq=(eq.lx+eq.lq)*(1-eqc2.p_lx);

eqc2.le=eq.le;


% Arrival rate
eqc2.x=par.chi*eqc2.lx^par.alpha_x;

%Quality
eqc2.dq=par.lambda*eqc2.lq^par.alpha_q;

%Entry
eqc2.xe=par.chi_e*eqc2.le^par.alpha_e;
eqc2.dqe=par.lambda_e;

%Growth
eqc2.g=eqc2.xe*eqc2.dqe+eqc2.x*eqc2.dq;

%3. Changing proportion of speed relative to total labor (keeping le constant)

eqc3.p_lx=eq0.lx/par0.L_I;

eqc3.lx=par.L_I*eqc3.p_lx;
eqc3.lq=par.L_I-eq.le-eqc3.lx;

eqc3.le=eq.le;


% Arrival rate
eqc3.x=par.chi*eqc3.lx^par.alpha_x;

%Quality
eqc3.dq=par.lambda*eqc3.lq^par.alpha_q;

%Entry
eqc3.xe=par.chi_e*eqc3.le^par.alpha_e;
eqc3.dqe=par.lambda_e;

%Growth
eqc3.g=eqc3.xe*eqc3.dqe+eqc3.x*eqc3.dq;


%Table

input.tableCaption='Labor Reallocation';
input.tableLabel='cf_lab';

resize=1;

if resize==1
    resize_in='\resizebox{\textwidth}{!}{';
    resize_out='}';
else
    resize_in='';
    resize_out='';
end

% make table header lines:

header = '\begin{tabular}{lccccc}';

latex = {'\begin{table}[h!]';'\centering'; ['\caption{',input.tableCaption,'}','\label{table:',input.tableLabel,'}']; ...
    header ;'\hline \hline'};

% latex = {'\begin{table}[h!]';'\centering'; ...
%     resize_in; header ;'\hline \hline'};

%latex = {'\begin{table}[h!]';'\centering'; header ;'\hline \hline'};

%Content of table
latex=[latex;...
    ['\textbf{Variable} & \textbf{' num2str(par.byear) '}  & \textbf{' num2str(par.cyear) '} & \textbf{$\Delta$ Quantity} & \textbf{$\Delta \frac{l_x}{l_q} $} \\'] ; ...
    ['Speed (Incumbents) & '  num2str(round(eq.x,5)) '& ' num2str(round(eq0.x,5)) '& ' num2str(round(eqc1.x,5))   '& ' num2str(round(eqc3.x,5))   ' \\' ];...
    ['Quality (Incumbents) & '  num2str(round(eq.dq,5)) '& ' num2str(round(eq0.dq,5)) '& ' num2str(round(eqc1.dq,5)) '& ' num2str(round(eqc3.dq,5))    ' \\' ];...
    ['Speed Entrants & '  num2str(round(eq.xe,5)) '& ' num2str(round(eq0.xe,5)) '& ' num2str(round(eqc1.xe,5))  '& ' num2str(round(eqc3.xe,5))  ' \\' ];...
    ['Quality Entrants & '  num2str(round(eq.dqe,5)) '& ' num2str(round(eq0.dqe,5)) '& ' num2str(round(eqc1.dqe,5)) '& ' num2str(round(eqc3.dqe,5))  ' \\' ];...
    ['Growth & '  num2str(100*round(eq.g,6)) '\% & ' num2str(100*round(eq0.g,6)) '\% & ' num2str(100*round(eqc1.g,6)) '\% & ' num2str(100*round(eqc3.g,6))  '\% \\' ];...
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


%% Changing labor for growth
ng=100;

lemax=0.5*par.L_I;
lemax0=0.5*par0.L_I;

%------------------------------
%1. Reallocating entrant labor
%------------------------------

leg=linspace(0,lemax,ng);
lincg=par.L_I-leg;

le0g=linspace(0,lemax0,ng);
linc0g=par0.L_I-le0g;

pleg=leg/par.L_I;
ple0g=le0g/par0.L_I;

plx=eq.lx/(eq.lq+eq.lx);
plx0=eq0.lx/(eq0.lq+eq0.lx);

lxg=plx*lincg; lqg=(1-plx)*lincg;
lx0g=plx0*linc0g; lq0g=(1-plx0)*linc0g;

%Arrival rate
xg=par.chi*lxg.^par.alpha_x;
x0g=par0.chi*lx0g.^par0.alpha_x;

%Quality
dqg=par.lambda*lqg.^par.alpha_q;
dq0g=par0.lambda*lq0g.^par0.alpha_q;

%Entry
xeg=par.chi_e*leg.^par.alpha_e;
dqeg=par.lambda_e;

xe0g=par0.chi_e*le0g.^par0.alpha_e;
dqe0g=par0.lambda_e;

%Growth
gg=xeg.*dqeg+xg.*dqg;
g0g=xe0g.*dqe0g+x0g.*dq0g;


%Plots

%For latex interpreter
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

%Speed
fig.xcle=figure; hold all
pb=plot(pleg,xg); pb.Color=par.opt.light_blue;
p0=plot(ple0g,x0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.le/par.L_I,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.le/par0.L_I,eq0.x); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Entrants')
ylabel('Arrival Rate')

%Quality
fig.dqcle=figure; hold all
pb=plot(pleg,dqg); pb.Color=par.opt.light_blue;
p0=plot(pleg,dq0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.le/par.L_I,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.le/par0.L_I,eq0.dq); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Entrants')
ylabel('Quality')

%Growth
fig.gcle=figure; hold all
pb=plot(pleg,100*gg); pb.Color=par.opt.light_blue;
p0=plot(pleg,100*g0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.le/par.L_I,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.le/par0.L_I,100*eq0.g); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Entrants')
ylabel('Growth \%')


%-------------------------------------------
%2. Changing the speed and quality of labor
%-------------------------------------------

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

%From estimation
x=par.chi*eq.lx.^par.alpha_x;
x0=par0.chi*eq0.lx.^par0.alpha_x;

dq=par.lambda*eq.lq.^par.alpha_q;
dq0=par0.lambda*eq0.lq.^par0.alpha_q;

xe=par.chi_e*eq.le.^par.alpha_e;
xe0=par0.chi_e*eq0.le.^par0.alpha_e;

dqe=par.lambda_e;
dqe0=par0.lambda_e;

g=xe*dqe+x*dq;
g0=xe0*dqe0+x0*dq0;

%Speed
fig.xclxq=figure; hold all
pb=plot(plxg,xg); pb.Color=par.opt.light_blue;
p0=plot(plxg,x0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.lx/linc,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.lx/linc0,eq0.x); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion To Speed')
ylabel('Arrival Rate')

%Quality
fig.dqclxq=figure; hold all
pb=plot(plxg,dqg); pb.Color=par.opt.light_blue;
p0=plot(plxg,dq0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.lx/linc,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.lx/linc0,eq0.dq); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion To Speed')
ylabel('Quality')

%Growth
fig.gclxq=figure; hold all
pb=plot(plxg,100*gg); pb.Color=par.opt.light_blue;
p0=plot(plxg,100*g0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.lx/linc,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.lx/linc0,100*eq0.g); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
%ylim([0.5,1.8])
ylim([1,1.8])
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion To Speed')
ylabel('Growth \%')


%-------------------------------------------
%2b. Changing the speed and quality of labor (Relative to total labor)
%-------------------------------------------

ng=1000;

%Incumbent labor
linc=eq.lx+eq.lq;
linc0=eq0.lx+eq0.lq;

% lxg=linspace(0.01*par.L_I,0.99*par.L_I,ng);
% lx0g=linspace(0.01*par0.L_I,0.99*par0.L_I,ng);

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

%From estimation
x=par.chi*eq.lx.^par.alpha_x;
x0=par0.chi*eq0.lx.^par0.alpha_x;

dq=par.lambda*eq.lq.^par.alpha_q;
dq0=par0.lambda*eq0.lq.^par0.alpha_q;

xe=par.chi_e*eq.le.^par.alpha_e;
xe0=par0.chi_e*eq0.le.^par0.alpha_e;

dqe=par.lambda_e;
dqe0=par0.lambda_e;

g=xe*dqe+x*dq;
g0=xe0*dqe0+x0*dq0;

%Proportion relative to total

%Speed
fig.xclxq_tot=figure; hold all
pb=plot(lxg/par.L_I,xg); pb.Color=par.opt.light_blue;
p0=plot(lx0g/par0.L_I,x0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.lx/par.L_I,eq.x); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.lx/par0.L_I,eq0.x); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion To Speed')
ylabel('Arrival Rate')

%Quality
fig.dqclxq_tot=figure; hold all
pb=plot(lxg/par.L_I,dqg); pb.Color=par.opt.light_blue;
p0=plot(lx0g/par0.L_I,dq0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.lx/par.L_I,eq.dq); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.lx/par0.L_I,eq0.dq); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion To Speed')
ylabel('Quality')

%Growth
fig.gclxq_tot=figure; hold all
pb=plot(lxg/par.L_I,100*gg); pb.Color=par.opt.light_blue;
p0=plot(lx0g/par0.L_I,100*g0g,'--'); p0.Color=par.opt.maroon;
sb=scatter(eq.lx/par.L_I,100*eq.g); sb.Marker='o'; sb.MarkerFaceColor='w'; sb.MarkerEdgeColor=par.opt.light_blue;
s0=scatter(eq0.lx/par0.L_I,100*eq0.g); s0.Marker='D'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.maroon;
%ylim([0.5,1.8])
ylim([1,1.8])
lgd=legend([pb,p0],{'2010 Parameters','1980 Parameters'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion To Speed')
ylabel('Growth \%')

%-----------------------------
%Optimal growth (Counterfactual)
%-----------------------------

%Find max
[~,indmax]=max(gg);

%Speed
fig.xcmaxlxq=figure; hold all
p0=plot(plxg,100*xg); p0.Color=par.opt.light_blue;
s0=scatter(eq.lx/linc,100*eq.x); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.light_blue;
smax=scatter(plxg(indmax),100*xg(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=par.opt.maroon; smax.MarkerEdgeColor=par.opt.maroon;
lgd=legend([s0,smax],{['Calibration ' num2str(par.iyear)],'Optimal'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Labor To Speed')
ylabel('Arrival Rate')


%Quality
fig.dqcmaxlxq=figure; hold all
p0=plot(plxg,100*dqg); p0.Color=par.opt.light_blue;
s0=scatter(eq.lx/linc,100*eq.dq); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.light_blue;
smax=scatter(plxg(indmax),100*dqg(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=par.opt.maroon; smax.MarkerEdgeColor=par.opt.maroon;
lgd=legend([s0,smax],{['Calibration ' num2str(par.iyear)],'Optimal'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Labor To Speed')
ylabel('Quality')


%Growth
fig.gcmaxlxq=figure; hold all
p0=plot(plxg,100*gg); p0.Color=par.opt.light_blue;
s0=scatter(eq.lx/linc,100*eq.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.light_blue;
smax=scatter(plxg(indmax),100*gg(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=par.opt.maroon; smax.MarkerEdgeColor=par.opt.maroon;
ylim([1.1,1.4])
lgd=legend([s0,smax],{['Calibration ' num2str(par.iyear)],'Optimal'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Labor To Speed')
ylabel('Growth \%')

%-----------------------------
%Optimal growth (Counterfactual) (relative to total)
%-----------------------------

%Find max
[~,indmax]=max(gg);

%Speed
fig.xcmaxlxq_tot=figure; hold all
p0=plot(lxg/par.L_I,100*xg); p0.Color=par.opt.light_blue;
s0=scatter(eq.lx/par.L_I,100*eq.x); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.light_blue;
smax=scatter(lxg(indmax)/par.L_I,100*xg(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=par.opt.maroon; smax.MarkerEdgeColor=par.opt.maroon;
lgd=legend([s0,smax],{['Calibration ' num2str(par.iyear)],'Optimal'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Labor To Speed')
ylabel('Arrival Rate')


%Quality
fig.dqcmaxlxq_tot=figure; hold all
p0=plot(lxg/par.L_I,100*dqg); p0.Color=par.opt.light_blue;
s0=scatter(eq.lx/par.L_I,100*eq.dq); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.light_blue;
smax=scatter(lxg(indmax)/par.L_I,100*dqg(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=par.opt.maroon; smax.MarkerEdgeColor=par.opt.maroon;
lgd=legend([s0,smax],{['Calibration ' num2str(par.iyear)],'Optimal'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Labor To Speed')
ylabel('Quality')


%Growth
fig.gcmaxlxq_tot=figure; hold all
p0=plot(lxg/par.L_I,100*gg); p0.Color=par.opt.light_blue;
s0=scatter(eq.lx/par.L_I,100*eq.g); s0.Marker='o'; s0.MarkerFaceColor='w'; s0.MarkerEdgeColor=par.opt.light_blue;
smax=scatter(lxg(indmax)/par.L_I,100*gg(indmax),80);  smax.Marker='pentagram'; smax.MarkerFaceColor=par.opt.maroon; smax.MarkerEdgeColor=par.opt.maroon;
ylim([1.1,1.4])
lgd=legend([s0,smax],{['Calibration ' num2str(par.iyear)],'Optimal'});lgd.EdgeColor="none"; lgd.Location='best';
xlabel('Proportion of Labor To Speed')
ylabel('Growth \%')



%Save figures
%par.opt.save_fig=1;
if ~isfield(par,'save_str')
    par.save_str='';
end

if par.opt.save_fig==1
    exportgraphics(fig.xcle,[par.opt.dir_fig 'xcle_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.dqcle,[par.opt.dir_fig 'dqcle_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.gcle,[par.opt.dir_fig 'gcle_intq_B' par.save_str '.pdf'])

    exportgraphics(fig.xclxq,[par.opt.dir_fig 'xclxq_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.dqclxq,[par.opt.dir_fig 'dqclxq_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.gclxq,[par.opt.dir_fig 'gclxq_intq_B' par.save_str '.pdf'])

    exportgraphics(fig.xclxq_tot,[par.opt.dir_fig 'xclxq_tot_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.dqclxq_tot,[par.opt.dir_fig 'dqclxq_tot_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.gclxq_tot,[par.opt.dir_fig 'gclxq_tot_intq_B' par.save_str '.pdf'])


    %Extra counterfactual
    exportgraphics(fig.xcmaxlxq,[par.opt.dir_fig 'xcmaxlxq_intq_B'  par.save_str '.pdf'])
    exportgraphics(fig.dqcmaxlxq,[par.opt.dir_fig 'dqcmaxlxq_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.gcmaxlxq,[par.opt.dir_fig 'gcmaxlxq_intq_B' par.save_str '.pdf'])

    exportgraphics(fig.xcmaxlxq_tot,[par.opt.dir_fig 'xcmaxlxq_tot_intq_B'  par.save_str '.pdf'])
    exportgraphics(fig.dqcmaxlxq_tot,[par.opt.dir_fig 'dqcmaxlxq_tot_intq_B' par.save_str '.pdf'])
    exportgraphics(fig.gcmaxlxq_tot,[par.opt.dir_fig 'gcmaxlxq_tot_intq_B' par.save_str '.pdf'])
end

 

end