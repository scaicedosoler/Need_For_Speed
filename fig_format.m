%Function for the format of the figures

function par=fig_format()

%Formatting figures
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
% set(0,'DefaultTextFontname', 'Garamond')
% set(0,'DefaultAxesFontName', 'Garamond')
set(0,'defaultAxesFontSize', 15)
set(0,'defaultlinelinewidth', 1.5)
set(0,'defaulttextInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex')

%Plotting options
par.opt.light_blue=[0 0.4470 0.7410];
par.opt.maroon=[0.6350, 0.0780, 0.1840];
par.opt.green=[0.4660, 0.6740, 0.1880];
par.opt.orange=[255,128,0]/255;
par.opt.lime=[0,255,0]/255;


end