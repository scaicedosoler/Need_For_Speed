%Function that fits moments for no endogenous quality model

function smm=fit_fun_mom_intq_B_nq(par,dmom)

%Choose moments to match: mom=(g)
par.nmom=1;
par.wmom=ones(1,par.nmom);

lb=0;
ub=10;

%Initial parameters
par0v(1)=par.lambda;

scoref=@(parv) score_mom_intq_B_nq(parv,dmom,par)*par.wmom';
scoref(par0v)

%Optimization
[smm.parv_sol, smm.score] =  patternsearch(scoref,par0v,[],[],[],[],lb,ub);

%Solution in structure
[smm.f,smm.par,smm.eq,smm.mom]=score_mom_intq_B_nq(smm.parv_sol,dmom,par);

end