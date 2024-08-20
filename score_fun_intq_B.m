%Score function given moment weights
function [score,par,eq]=score_fun_intq_B(parv,dmom,par)

if ~isfield(par,'wmom')
    npar=length(parv);
    par.wmom=ones(1,npar); %Default weights
end


%Compute score
[f,par,eq,mom]=score_mom_intq_B(parv,dmom,par);

score=f*par.wmom'

%Save minimum score

if score<0
    save([par.opt.dir_mat 'par_score_err.mat'],'par')
    error('Negative score value');
end

if par.opt.save_score==1

    if ~isfield(par.opt,'save_str')
        par.opt.save_str=[];
    end
    
    try
    load([ par.opt.dir_mat 'min_score_' par.opt.save_str '.mat'],'min_score')

    catch
        min_score=100;
    end

    if score<min_score
        min_score=score;
        save([ par.opt.dir_mat 'min_score_' par.opt.save_str '.mat'],'min_score','par','eq','mom')
    end

end

if par.opt.save_last_run==1
    %Save also latest run
    save([ par.opt.dir_mat 'last_run_' par.opt.save_str '.mat'],'score','par','eq')
end

end