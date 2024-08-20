%Counterfactual: Optimal Growth Comparison
function og=cf_opt_g_intq_B(par)


%Number of points for optimal growth
og.np=9;

%chi
    disp('-------------------------------------------------')
    disp('Optimal Growth: chi')
    disp('-------------------------------------------------')
    og.par='chi'; og.dsub=0;
    og_min=0.1; og_max=2; step=(og_max-og_min)/(og.np-1);
    og.chig=unique(sort([og_min:step:og_max,par.chi]));
    og.chi=og_fun_intq_B(par,og);

%alpha_x
    disp('-------------------------------------------------')
    disp('Optimal Growth: alpha_x')
    disp('-------------------------------------------------')
    og.par='alpha_x'; og.dsub=1;
    og_min=0; og_max=1-par.alpha_q-0.01; step=(og_max-og_min)/(og.np-1);
    og.alpha_xg=unique(sort([og_min:step:og_max,par.alpha_x]));
    og.alpha_x=og_fun_intq_B(par,og);

%alpha_q
    disp('-------------------------------------------------')
    disp('Optimal Growth: alpha_q')
    disp('-------------------------------------------------')
    og.par='alpha_q'; og.dsub=1;
    og_min=0; og_max=1-par.alpha_x-0.01; step=(og_max-og_min)/(og.np-1);
    og.alpha_qg=unique(sort([og_min:step:og_max,par.alpha_q]));
    og.alpha_q=og_fun_intq_B(par,og);



%lambda
    disp('-------------------------------------------------')
    disp('Optimal Growth: lambda')
    disp('-------------------------------------------------')
    og.par='lambda'; og.dsub=0;
    og_min=0.5; og_max=3; step=(og_max-og_min)/(og.np-1);
    og.lambdag=unique(sort([og_min:step:og_max,par.lambda]));
    og.lambda=og_fun_intq_B(par,og);

%chi_e
    disp('-------------------------------------------------')
    disp('Optimal Growth: chi_e')
    disp('-------------------------------------------------')
    og.par='chi_e'; og.dsub=1;
    og_min=0.001; og_max=0.1; step=(og_max-og_min)/(og.np-1);
    og.chi_eg=unique(sort([og_min:step:og_max,par.chi_e]));
    og.chi_e=og_fun_intq_B(par,og);

%lambda_e
    disp('-------------------------------------------------')
    disp('Optimal Growth: lambda_e')
    disp('-------------------------------------------------')
    og.par='lambda_e'; og.dsub=1;
    og_min=0.5; og_max=2; step=(og_max-og_min)/(og.np-1);
    og.lambda_eg=unique(sort([og_min:step:og_max,par.lambda_e]));
    og.lambda_e=og_fun_intq_B(par,og);

%alpha_e
    disp('-------------------------------------------------')
    disp('Optimal Growth: alpha_e')
    disp('-------------------------------------------------')
    og.par='alpha_e'; og.dsub=1;
    og_min=0.1; og_max=0.8; step=(og_max-og_min)/(og.np-1);
    og.alpha_eg=unique(sort([og_min:step:og_max,par.alpha_e]));
    og.alpha_e=og_fun_intq_B(par,og);


%chi_b
    disp('-------------------------------------------------')
    disp('Optimal Growth: chi_b')
    disp('-------------------------------------------------')
    og.par='chi_b'; og.dsub=1;
    og_min=0; og_max=40; step=(og_max-og_min)/(og.np-1);
    og.chi_bg=unique(sort([og_min:step:og_max,par.chi_b]));
    og.chi_b=og_fun_intq_B(par,og);

%gamma_q
    disp('-------------------------------------------------')
    disp('Optimal Growth: gamma_q')
    disp('-------------------------------------------------')
    og.par='gamma_q'; og.dsub=1;
    og_min=-0.2; og_max=0.2; step=(og_max-og_min)/(og.np-1);
    og.gamma_qg=unique(sort([og_min:step:og_max,par.gamma_q]));
    og.gamma_q=og_fun_intq_B(par,og);

%nu
    disp('-------------------------------------------------')
    disp('Optimal Growth: nu')
    disp('-------------------------------------------------')
    og.par='nu'; og.dsub=0;
    og_min=0.1; og_max=0.9; step=(og_max-og_min)/(og.np-1);
    og.nug=unique(sort([og_min:step:og_max,par.nu]));
    og.nu=og_fun_intq_B(par,og);

end